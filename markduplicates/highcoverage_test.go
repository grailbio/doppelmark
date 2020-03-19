package markduplicates

import (
	"fmt"
	"path/filepath"
	"testing"

	gbam "github.com/grailbio/bio/encoding/bam"
	"github.com/grailbio/bio/encoding/bamprovider"
	"github.com/grailbio/hts/sam"
	"github.com/grailbio/testutil"
	"github.com/stretchr/testify/assert"
)

func TestHighCoverage(t *testing.T) {
	ref1, _ := sam.NewReference("ref1", "", "", 3, nil, nil)
	ref2, _ := sam.NewReference("ref2", "", "", 3, nil, nil)
	header, _ := sam.NewHeader(nil, []*sam.Reference{ref1, ref2})
	assert.NotNil(t, header)

	shard0 := gbam.Shard{
		StartRef: ref1,
		EndRef:   ref1,
		Start:    0,
		End:      2,
		StartSeq: 0,
		EndSeq:   0,
		Padding:  0,
		ShardIdx: 0,
	}
	shard1 := gbam.Shard{
		StartRef: ref1,
		EndRef:   ref2,
		Start:    2,
		End:      1,
		StartSeq: 0,
		EndSeq:   0,
		Padding:  0,
		ShardIdx: 1,
	}
	shard2 := gbam.Shard{
		StartRef: ref2,
		EndRef:   ref2,
		Start:    1,
		End:      3,
		StartSeq: 0,
		EndSeq:   0,
		Padding:  0,
		ShardIdx: 2,
	}

	testCases := []struct {
		name                     string
		shard                    gbam.Shard
		records                  []*sam.Record
		expectedCoverageCounts   map[int][]int
		expectedHighCovIntervals []coverageInterval
	}{
		{
			name:  "shard0-only",
			shard: shard0,
			records: []*sam.Record{
				NewRecord("A", ref1, 0, r1F, 10, ref1, cigar2M),
			},
			expectedCoverageCounts: map[int][]int{
				0: []int{1, 1, 0},
				1: []int{0, 0, 0},
			},
			expectedHighCovIntervals: []coverageInterval{},
		},
		{
			name:  "shard0-partial",
			shard: shard0,
			records: []*sam.Record{
				NewRecord("A", ref1, 1, r1F, 10, ref1, cigar2M),
			},
			expectedCoverageCounts: map[int][]int{
				0: []int{0, 1, 0},
				1: []int{0, 0, 0},
			},
			expectedHighCovIntervals: []coverageInterval{},
		},
		{
			name:  "shard1-partial",
			shard: shard1,
			records: []*sam.Record{
				NewRecord("A", ref1, 2, r1F, 10, ref1, cigar2M),
			},
			expectedCoverageCounts: map[int][]int{
				0: []int{0, 0, 1},
				1: []int{0, 0, 0},
			},
			expectedHighCovIntervals: []coverageInterval{},
		},
		{
			name:  "shard1-partial2",
			shard: shard1,
			records: []*sam.Record{
				NewRecord("A", ref2, 0, r1F, 10, ref1, cigar2M),
			},
			expectedCoverageCounts: map[int][]int{
				0: []int{0, 0, 0},
				1: []int{1, 0, 0},
			},
			expectedHighCovIntervals: []coverageInterval{},
		},
		{
			name:  "shard2-starts-before-shard",
			shard: shard2,
			records: []*sam.Record{
				NewRecord("A", ref2, 0, r1F, 10, ref1, cigar2M),
			},
			expectedCoverageCounts: map[int][]int{
				0: []int{0, 0, 0},
				1: []int{0, 1, 0},
			},
			expectedHighCovIntervals: []coverageInterval{},
		},
		{
			name:  "shard2-inshard",
			shard: shard2,
			records: []*sam.Record{
				NewRecord("A", ref2, 1, r1F, 10, ref1, cigar2M),
			},
			expectedCoverageCounts: map[int][]int{
				0: []int{0, 0, 0},
				1: []int{0, 1, 1},
			},
			expectedHighCovIntervals: []coverageInterval{},
		},
		{
			name:  "shard2-partial",
			shard: shard2,
			records: []*sam.Record{
				NewRecord("A", ref2, 2, r1F, 10, ref1, cigar2M),
			},
			expectedCoverageCounts: map[int][]int{
				0: []int{0, 0, 0},
				1: []int{0, 0, 1},
			},
			expectedHighCovIntervals: []coverageInterval{},
		},
		{
			name:  "shard0-two",
			shard: shard0,
			records: []*sam.Record{
				NewRecord("A", ref1, 0, r1F, 10, ref1, cigar2M),
				NewRecord("A", ref1, 1, r1F, 10, ref1, cigar2M),
			},
			expectedCoverageCounts: map[int][]int{
				0: []int{1, 2, 0},
				1: []int{0, 0, 0},
			},
			expectedHighCovIntervals: []coverageInterval{
				coverageInterval{
					refId:        0,
					start:        1,
					end:          2,
					meanCoverage: 2.0,
				},
			},
		},
		{
			name:  "shard1-two",
			shard: shard1,
			records: []*sam.Record{
				NewRecord("A", ref1, 1, r1F, 10, ref1, cigar2M),
				NewRecord("A", ref1, 2, r1F, 10, ref1, cigar2M),
			},
			expectedCoverageCounts: map[int][]int{
				0: []int{0, 0, 2},
				1: []int{0, 0, 0},
			},
			expectedHighCovIntervals: []coverageInterval{
				coverageInterval{
					refId:        0,
					start:        2,
					end:          3,
					meanCoverage: 2.0,
				},
			},
		},
		{
			name:  "shard2-two",
			shard: shard2,
			records: []*sam.Record{
				NewRecord("A", ref2, 1, r1F, 10, ref1, cigar2M),
				NewRecord("A", ref2, 2, r1F, 10, ref1, cigar2M),
			},
			expectedCoverageCounts: map[int][]int{
				0: []int{0, 0, 0},
				1: []int{0, 1, 2},
			},
			expectedHighCovIntervals: []coverageInterval{
				coverageInterval{
					refId:        1,
					start:        2,
					end:          3,
					meanCoverage: 2.0,
				},
			},
		},
	}

	for _, testCase := range testCases {
		t.Run(testCase.name, func(t *testing.T) {
			// References ref1 and ref2
			coverageCounts := map[int][]int{
				0: make([]int, ref1.Len()),
				1: make([]int, ref2.Len()),
			}
			c := coverageCalculator{
				coverageCounts: &coverageCounts,
			}
			for _, r := range testCase.records {
				err := c.Process(testCase.shard, r)
				assert.NoError(t, err)
			}
			assert.Equal(t, testCase.expectedCoverageCounts, coverageCounts)

			// identify high-coverage intervals
			highCovIntervals := getHighCoverageIntervals(coverageCounts, 1)
			assert.Equal(t, testCase.expectedHighCovIntervals, highCovIntervals)
		})
	}
}

func TestGetHighCoverageIntervals(t *testing.T) {
	testCases := []struct {
		name        string
		coverage    map[int][]int
		maxCoverage int
		expected    []coverageInterval
	}{
		{
			name: "basic",
			coverage: map[int][]int{
				0: []int{0, 0, 1, 2, 3},
				1: []int{2, 2, 1, 3},
				2: []int{1, 1, 4, 2, 1},
				3: []int{1, 1, 4, 1, 1},
			},
			maxCoverage: 1,
			expected: []coverageInterval{
				coverageInterval{
					refId:        0,
					start:        3,
					end:          5,
					meanCoverage: 2.5,
				},
				coverageInterval{
					refId:        1,
					start:        0,
					end:          2,
					meanCoverage: 2,
				},
				coverageInterval{
					refId:        1,
					start:        3,
					end:          4,
					meanCoverage: 3,
				},
				coverageInterval{
					refId:        2,
					start:        2,
					end:          4,
					meanCoverage: 3,
				},
				coverageInterval{
					refId:        3,
					start:        2,
					end:          3,
					meanCoverage: 4,
				},
			},
		},
	}

	for _, testCase := range testCases {
		t.Run(testCase.name, func(t *testing.T) {
			highCovIntervals := getHighCoverageIntervals(testCase.coverage, testCase.maxCoverage)
			assert.Equal(t, testCase.expected, highCovIntervals)
		})
	}
}

func TestIsInHighCoverageShard(t *testing.T) {
	highCovMap := getCoverageMap([]coverageInterval{
		coverageInterval{
			refId:        0,
			start:        22,
			end:          23,
			meanCoverage: 5,
		},
		coverageInterval{
			refId:        1,
			start:        43,
			end:          45,
			meanCoverage: 10,
		},
	})

	tests := []struct {
		record    *sam.Record
		assertion assert.BoolAssertionFunc
	}{
		{
			NewRecord("A:::1:10:1:1", chr1, 0, r1F, 22, chr1, cigar0),
			assert.True,
		}, {
			NewRecord("A:::1:10:1:1", chr1, 0, r1F, 44, chr2, cigar0),
			assert.True,
		}, {
			NewRecord("A:::1:10:1:1", chr1, 22, r1F, 22, chr1, cigar0),
			assert.True,
		}, {
			NewRecord("A:::1:10:1:1", chr1, 22, r1F, 44, chr2, cigar0),
			assert.True,
		}, {
			NewRecord("A:::1:10:1:1", chr1, 22, r1F, 0, chr2, cigar0),
			assert.True,
		}, {
			NewRecord("A:::1:10:1:1", chr1, 22, r1F, 100, chr1, cigar0),
			assert.True,
		}, {
			NewRecord("A:::1:10:1:1", chr1, 90, r1F, 100, chr1, cigar0),
			assert.False,
		},
	}
	for _, test := range tests {
		found, _ := recOrMateInHighCovInterval(highCovMap, test.record)
		test.assertion(t, found)
	}
}

// Test end-to-end high-coverage removal using Mark(). This test
// creates two high coverage regions and then checks that the output
// has those two regions subsampled down to approximately the right
// number of reads.
func TestSubsampleCoverageMax(t *testing.T) {
	tempDir, cleanup := testutil.TempDir(t, "", "")
	defer cleanup()

	const (
		numRecords  = 10000
		coverageMax = 3000 // Subsample 0.1 of reads C and D.
	)

	outputPath := filepath.Join(tempDir, "foo.bam")
	opts := Opts{
		ShardSize:            100,
		Padding:              10,
		Parallelism:          1,
		QueueLength:          10,
		EmitUnmodifiedFields: true,
		Format:               "bam",
		OutputPath:           outputPath,
		CoverageMax:          coverageMax,
		Seed:                 1233,
	}

	var records []*sam.Record
	records = append(records, NewRecordSeq("A", chr1, 5, r1F, 5, chr1, cigar2M, "AC", "FF"))
	records = append(records, NewRecordSeq("A", chr1, 5, r2R, 5, chr1, cigar2M, "AC", "FF"))
	records = append(records, NewRecordSeq("B", chr1, 10, r1F, 10, chr1, cigar2M, "AC", "FF"))
	records = append(records, NewRecordSeq("B", chr1, 10, r2R, 10, chr1, cigar2M, "AC", "FF"))

	// B, C_i, and D_i overlap and create a region of meanCoverage=30001 at chr1:11-13
	for i := 0; i < numRecords; i++ {
		records = append(records, NewRecordSeq(fmt.Sprintf("C%d", i), chr1, 11, r1F, 11, chr1, cigar2M, "AC", "FF"))
		records = append(records, NewRecordSeq(fmt.Sprintf("C%d", i), chr1, 11, r2R, 11, chr1, cigar2M, "AC", "FF"))
		records = append(records, NewRecordSeq(fmt.Sprintf("D%d", i), chr1, 11, r1F, 100, chr2, cigar2M, "AC", "FF"))
	}
	records = append(records, NewRecordSeq("E", chr1, 15, r1F, 15, chr1, cigar2M, "AC", "FF"))
	records = append(records, NewRecordSeq("E", chr1, 15, r2R, 15, chr1, cigar2M, "AC", "FF"))
	// The R2 for D_i also creates a high-coverage region at chr1:100-103 but is smaller than the one at chr1:11-13.
	for i := 0; i < numRecords; i++ {
		records = append(records, NewRecordSeq(fmt.Sprintf("D%d", i), chr2, 100, r2R, 11, chr1, cigar2M, "AC", "FF"))
	}
	provider := bamprovider.NewFakeProvider(header, records)

	markDuplicates := &MarkDuplicates{
		Provider: provider,
		Opts:     &opts,
	}
	_, err := markDuplicates.Mark(nil)
	assert.NoError(t, err)
	for i, r := range records {
		t.Logf("input[%v]: %v begin %d end %d", i, r, r.Start(), r.End())
	}
	actualRecords := ReadRecords(t, outputPath)

	counts := make(map[string]int)
	for _, r := range actualRecords {
		counts[r.Name[0:1]]++
	}

	// counts for A, B, and E should be exact and not subsampled.
	assert.Equal(t, counts["A"], 2)
	assert.Equal(t, counts["B"], 2)
	assert.Equal(t, counts["E"], 2)

	// Check that the counts for C and D are within 10% of expected value
	expectedCount := float64(2*numRecords) * (float64(coverageMax) / (3 * numRecords))
	assert.Greater(t, float64(counts["C"]), expectedCount*0.9)
	assert.Less(t, float64(counts["C"]), expectedCount*1.1)
	assert.Greater(t, float64(counts["D"]), expectedCount*0.9)
	assert.Less(t, float64(counts["D"]), expectedCount*1.1)
}
