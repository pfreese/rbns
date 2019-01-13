package rbns

import (
	"errors"
	"gotest.tools/assert"
	"reflect"
	"testing"
)

func TestCountKmers(t *testing.T) {
	tables := []struct {
		s string
		k int
		kmerCounts KmerCountMap
	}{
		{"ABCDE", 2, KmerCountMap{"AB": 1, "BC": 1, "CD": 1, "DE": 1}},
	}

	for _, table := range tables {
		actual := CountKmers(table.s, table.k)
		if !reflect.DeepEqual(actual, table.kmerCounts) {
			t.Errorf("%dmer counts of %s was incorrect, got: %v, want: %v.",
				table.k, table.s, actual, table.kmerCounts)
		}
	}
}

func TestIsACGT(t *testing.T) {
	tables := []struct {
		s string
		exp bool // expected result
	}{
		{"A", true},
		{"C", true},
		{"G", true},
		{"T", true},
		{"U", false}, // no RNA
		{"", false},
		{"AA", false}, // only single characters expected
	}

	for _, table := range tables {
		actual := isACGT(table.s)
		if table.exp != actual {
			t.Errorf("Is %s one of A/C/G/T? was incorrect, got: %v, want: %v.",
				table.s, actual, table.exp)
		}
	}
}

func TestIsACGTOnlyKmer(t *testing.T) {
	tables := []struct {
		s string
		exp bool // expected result
	}{
		{"A", true},
		{"C", true},
		{"G", true},
		{"T", true},
		{"U", false}, // no RNA
		{"", false}, // must be at least length 1
		{"AA", true},
		{"TTAACAGATAACC", true},
		{" TTAACAGATAACC", false}, // extra space makes this invalid
		{"tTTAACAGATAACC", false}, // lower case t makes this invalid
	}

	for _, table := range tables {
		actual := isACGTonlyKmer(table.s)
		if table.exp != actual {
			t.Errorf("Is %s an A/C/G/T-only kmer? was incorrect, got: %v, want: %v.",
				table.s, actual, table.exp)
		}
	}
}

func TestFilterToACGTkmers(t *testing.T) {
	tables := []struct {
		inMap KmerCountMap
		filtMap KmerCountMap // expected result
	}{
		{KmerCountMap{"AC": 1, "BC": 1, "CD": 1, "DE": 1},
			KmerCountMap{"AC": 1},
		},
		{KmerCountMap{"GCCG": 10, "AAAA": 3, "CG": 0, "DE": 10},
			KmerCountMap{"GCCG": 10, "AAAA": 3, "CG": 0},
		},
		{KmerCountMap{},
			KmerCountMap{},
		},
	}

	for _, table := range tables {
		actual := FilterToACGTkmers(table.inMap)
		if !reflect.DeepEqual(actual, table.filtMap) {
			t.Errorf("Filtering of %v was incorrect, got: %v, want: %v.",
				table.inMap, actual, table.filtMap)
		}
	}
}

func TestCountsMapToFreqMap(t *testing.T) {
	tables := []struct {
		inMap KmerCountMap
		freqMap KmerFreqMap // expected result
		error error
	}{
		{KmerCountMap{"GCCG": 10, "AAAA": 5, "CG": 5},
			KmerFreqMap{"GCCG": 0.5, "AAAA": 0.25, "CG": 0.25},
			nil,
		},
		// Error if non-ACGT kmers
		{KmerCountMap{"GCCG": 10, "AAAA": 5, "CG": 5, "H": 10},
			nil,
			errors.New("KmerCountMap contains non-ACGT kmer H"),
		},
		// Error if negative kmer count
		{KmerCountMap{"GCCG": 10, "AAAA": 5, "CG": -1},
			nil,
			errors.New("KmerCountMap contains negative counts for"),
		},
		// Empty map returns empty map
		{KmerCountMap{},
			KmerFreqMap{},
			nil,
		},
	}

	for _, table := range tables {
		actual, err := CountsMapToFreqMap(table.inMap)
		// If the expected error is not nil,
		if table.error != nil {
			assert.ErrorContains(t, err, table.error.Error())
		} else {
			if !reflect.DeepEqual(actual, table.freqMap) {
				t.Errorf("Frequencies from %v was incorrect, got: %v, want: %v.",
					table.inMap, actual, table.freqMap)
			}
		}
	}
}