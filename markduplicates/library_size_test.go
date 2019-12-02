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
	"testing"

	"github.com/stretchr/testify/assert"
)

func TestEstimateLibrarySize(t *testing.T) {
	tests := []struct {
		readPairs       uint64
		uniqueReadPairs uint64
		expected        uint64
	}{
		{1000000, 800000, 2154184},
		{171512300, 171512299, 14708234445116054},
	}

	for _, test := range tests {
		v, err := estimateLibrarySize(test.readPairs, test.uniqueReadPairs)
		assert.NoError(t, err)
		assert.InEpsilon(t, test.expected, v, 0.0000000001)
	}
}
