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

	"github.com/grailbio/base/log"
)

type Orientation uint8

const (
	f  = iota // Forward (single fragment)
	r  = iota // Reverse (single fragment)
	ff = iota // Forward, Forward
	fr = iota // Forward, Reverse
	rf = iota // Reverse, Forward
	rr = iota // Reverse, Reverse
)

// duplicateKey is a unique key for each group of duplicates.  If both
// left and right are populated, the left most unclipped 5' position will
// reside in left.  If only one read is populated, it will reside in left,
// and .isSingle() returns true.
type duplicateKey struct {
	leftRefId   int
	leftPos     int
	rightRefId  int
	rightPos    int
	Orientation Orientation
	Strand      strand
}

func (k *duplicateKey) String() string {
	return fmt.Sprintf("(%d,%d,%d,%d,0x%x,%d)", k.leftRefId, k.leftPos,
		k.rightRefId, k.rightPos, k.Orientation, k.Strand)
}

func (k *duplicateKey) isSingle() bool {
	return k.Orientation == f || k.Orientation == r
}

func orientationByteSingle(reversed bool) Orientation {
	if reversed {
		return r
	}
	return f
}

func leftOrientation(o Orientation) Orientation {
	if o == f || o == r {
		log.Fatal("expected pair orientation, got single fragment orientation")
		return r
	} else if o == ff || o == fr {
		return f
	} else {
		return r
	}
}

func rightOrientation(o Orientation) Orientation {
	if o == f || o == r {
		log.Fatal("expected pair orientation, got single fragment orientation")
		return r
	} else if o == ff || o == rf {
		return f
	} else {
		return r
	}
}

func orientationBytePair(leftReversed, rightReversed bool) Orientation {
	if leftReversed {
		if rightReversed {
			return rr
		}
		return rf
	}
	if rightReversed {
		return fr
	}
	return ff
}
