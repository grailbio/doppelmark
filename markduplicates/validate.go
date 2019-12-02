// Copyright 2019 Grail Inc.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
package markduplicates

import (
	"fmt"

	"github.com/grailbio/bio/encoding/bamprovider"
)

func validate(opts *Opts) error {
	if opts.BamFile == "" {
		return fmt.Errorf("you must specify a bam file with --bam")
	}
	if opts.ShardSize <= 0 {
		return fmt.Errorf("shard-size must be non-zero")
	}
	if opts.Padding < 0 {
		return fmt.Errorf("padding must be non-negative")
	}
	if opts.Padding >= opts.ShardSize {
		return fmt.Errorf("padding must be less than shard-size")
	}
	if opts.MinBases <= 0 {
		return fmt.Errorf("min-bases should be positive")
	}
	if opts.IndexFile == "" {
		opts.IndexFile = opts.BamFile + ".bai"
	}
	if len(opts.UmiFile) > 0 && !opts.UseUmis {
		return fmt.Errorf("umi-file is set, but use-umis is false")
	}
	if opts.ScavengeUmis > -1 && !opts.UseUmis {
		return fmt.Errorf("scavenge-umis is set, but use-umis is false")
	}
	if opts.ScavengeUmis > -1 && opts.UmiFile == "" {
		return fmt.Errorf("scavenge-umis is set, but umi-file is empty")
	}
	if bamprovider.ParseFileType(opts.Format) == bamprovider.Unknown {
		return fmt.Errorf("unknown outputformat %s", opts.Format)
	}
	return nil
}
