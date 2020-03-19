package markduplicates

import (
	"github.com/grailbio/base/intervalmap"
	"github.com/grailbio/base/log"
	"github.com/grailbio/bio/encoding/bam"
	"github.com/grailbio/hts/sam"
)

type coverageInterval struct {
	refId        int
	start        int
	end          int
	meanCoverage float64
}

// coverageCalculator calculates the per-base coverage from within GetDistantMates.
// It writes the coverage counts to coverageCounts.
type coverageCalculator struct {
	coverageCounts *map[int][]int
}

func (m *coverageCalculator) Process(shard bam.Shard, r *sam.Record) error {
	// Count the number of bases that precede the shard.
	basesPreShard := 0
	for p := r.Start(); p < r.End(); p++ {
		if !shard.CoordInShard(0, bam.NewCoord(r.Ref, p, 0)) {
			basesPreShard++
		} else {
			break
		}
	}
	if basesPreShard >= r.Len() {
		return nil
	}

	// Count the number of bases that actually overlap the shard.
	pos := r.Start()
	basesInShard := r.Len()
	for p := r.End() - 1; p >= pos; p-- {
		if !shard.CoordInShard(0, bam.NewCoord(r.Ref, p, 0)) {
			basesInShard--
		} else {
			break
		}
	}

	// Increment coverage counters for bases that overlap the shard.
	// Unmapped reads do not contribute to coverage counts.
	counted := 0
	offset := 0
	for _, co := range r.Cigar {
		if co.Type().Consumes().Reference == 1 {
			for i := 0; i < co.Len() && counted < basesInShard && pos+offset < r.Ref.Len(); i++ {
				if offset >= basesPreShard {
					(*m.coverageCounts)[r.Ref.ID()][pos+offset]++
					counted++
				}
				offset++
			}
		}
	}
	return nil
}

func (m *coverageCalculator) Close(_ bam.Shard) {}

// getHighCoverageIntervals takes the coverageCounts computed by coverageCalculator
// and returns a slice of coverageIntervals where the coverage is higher than maxCoverage.
// The output is sorted by refId and then position.
func getHighCoverageIntervals(coverage map[int][]int, maxCoverage int) []coverageInterval {
	highCovIntervals := make([]coverageInterval, 0)
	for refId := 0; refId < len(coverage); refId++ {
		refCoverage := coverage[refId]
		var start, end, total int
		for pos := range refCoverage {
			if refCoverage[pos] > maxCoverage {
				log.Printf("highcoverage ref %d pos %d depth %d", refId, pos, refCoverage[pos])
				if pos == 0 || (pos > 0 && refCoverage[pos-1] <= maxCoverage) {
					start = pos
					total = 0
				}
				total += refCoverage[pos]
				if pos == len(refCoverage)-1 {
					end = pos + 1
					highCovIntervals = append(highCovIntervals, coverageInterval{
						refId:        refId,
						start:        start,
						end:          end,
						meanCoverage: float64(total) / float64(end-start),
					})
					log.Printf("highcoverage range: %d %d-%d depth %f", refId, start, end,
						float64(total)/float64(end-start))
				}
			}
			if refCoverage[pos] <= maxCoverage {
				if pos > 0 && refCoverage[pos-1] > maxCoverage {
					end = pos
					highCovIntervals = append(highCovIntervals, coverageInterval{
						refId:        refId,
						start:        start,
						end:          end,
						meanCoverage: float64(total) / float64(end-start),
					})
					log.Printf("highcoverage range: %d %d-%d depth %f", refId, start, end,
						float64(total)/float64(end-start))
				}
			}
		}
	}
	return highCovIntervals
}

// coverageMap associates each refId to an intervalmap containing
// high-coverage intervals.
type coverageMap map[int]*intervalmap.T

// getCoverageMap returns a coverageMap that allows efficient
// intersection calls, given a slice of coverageIntervals.
func getCoverageMap(intervals []coverageInterval) coverageMap {
	allEntries := make(map[int][]intervalmap.Entry)
	for _, interval := range intervals {
		allEntries[interval.refId] = append(
			allEntries[interval.refId],
			intervalmap.Entry{
				Interval: intervalmap.Interval{
					Start: int64(interval.start),
					Limit: int64(interval.end),
				},
				Data: interval.meanCoverage,
			})
	}

	returnMap := make(coverageMap)
	for refId, entries := range allEntries {
		returnMap[refId] = intervalmap.New(entries)
	}
	return returnMap
}
