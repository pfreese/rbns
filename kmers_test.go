package rbns

import (
	"errors"
	"github.com/stretchr/testify/assert"
	"reflect"
	"testing"
)

// Test for panics based on: https://stackoverflow.com/questions/31595791/how-to-test-panics
func TestAllKmers(t *testing.T) {
	tables := []struct {
		name string // name of the test
		k int
		kmers KmersSlice // expected returned kmers
		wantPanic bool // whether the test should trigger a panic (i.e., k < 1)
	}{
		{"k=0",
			0,
			nil,
			true,},
		{"k=1",
			1,
			KmersSlice{"A", "C", "G", "T"},
			false,
		},
		{"k=2",
			2,
			KmersSlice{
				"AA", "AC", "AG", "AT",
				"CA", "CC", "CG", "CT",
				"GA", "GC", "GG", "GT",
				"TA", "TC", "TG", "TT",
			},
			false,
		},
	}

	for _, tt := range tables {
		// Test if there's a panic
		t.Run(tt.name, func(t *testing.T) {
			defer func() {
				r := recover()
				if (r != nil) != tt.wantPanic {
					t.Errorf("AllKmers() recover = %v, wantPanic = %v", r, tt.wantPanic)
				}
			}()
			// If no panic, test kmers are as expected
			got := AllKmers(tt.k)
			if !reflect.DeepEqual(got, tt.kmers) {
				t.Errorf("AllKmers() = %v, want %v", got, tt.kmers)
			}
		})
	}
}

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
			assert.Error(t, err)
			assert.Contains(t, err.Error(), table.error.Error())
		} else {
			if !reflect.DeepEqual(actual, table.freqMap) {
				t.Errorf("Frequencies from %v was incorrect, got: %v, want: %v.",
					table.inMap, actual, table.freqMap)
			}
		}
	}
}

func TestCheckValidKmerFreqMap(t *testing.T) {
	tables := []struct {
		freqMap KmerFreqMap
		k int // size of kmers in freqMap
		error error  // expected error
	}{
		// Valid
		{KmerFreqMap{"A": 0.5, "C": 0.5},
			1,
			nil,
		},
		// Wrong length kmers
		{KmerFreqMap{"A": 0.5, "C": 0.5},
			2,
			errors.New("KmerFreqMap contains non-length 2 kmer:"),
		},
		// Negative entry
		{KmerFreqMap{"A": -0.2, "C": 0.5},
			1,
			errors.New("freq for A is <0"),
		},
		// Does not sum to 1
		{KmerFreqMap{"A": 0.4999, "C": 0.5},
			1,
			errors.New("total freq for all kmers is not 1"),
		},
		// non-ACGT kmer
		{KmerFreqMap{"A": 0.5, "N": 0.5},
			1,
			errors.New("non-ACGT kmer: N"),
		},
	}

	for _, table := range tables {
		err := checkValidKmerFreqMap(table.freqMap, table.k)
		// If the expected error is not nil,
		if table.error != nil {
			assert.Error(t, err)
			assert.Contains(t, err.Error(), table.error.Error())
		} else {
			assert.NoError(t, err)
		}
	}
}

func TestGetK(t *testing.T) {
	tables := []struct {
		kmersMap KmerFreqMap
		k int // expected result
	}{
		{KmerFreqMap{"AC": 0.25, "GG": 0.25, "CG": 0.25, "GA": 0.25},
			2,
		},
		{KmerFreqMap{"GCCG": 0.1, "AAAA": 0.9},
			4,
		},
		{KmerFreqMap{},
			0,
		},
	}

	for _, table := range tables {
		actual := getk(table.kmersMap)
		if actual != table.k {
			t.Errorf("k from %v was incorrect, got: %v, want: %v.",
				table.kmersMap, actual, table.k)
		}
	}
}

func TestKmerEnrichments(t *testing.T) {
	tables := []struct {
		pdMap KmerFreqMap
		inputMap KmerFreqMap
		rMap KmerRMap // expected enrichments
		error error  // expected error
	}{
		// "A" and "C" have both positive frequencies
		{KmerFreqMap{"A": 0.5, "C": 0.5},
			KmerFreqMap{"A": 0.25, "C": 0.75},
			KmerRMap{"A": 0.5/0.25, "C": 0.5/0.75, "G": 1, "T": 1},
			nil,
		},
		// Only "A" has frequency in both libraries
		{KmerFreqMap{"A": 0.5, "C": 0.5},
			KmerFreqMap{"A": 0.25, "G": 0.75},
			KmerRMap{"A": 0.5/0.25, "C": 1, "G": 1, "T": 1},
			nil,
		},
		// Test a couple of errors:
		// If the PD freqs do not sum to 1
		{KmerFreqMap{"A": 0.5, "C": 1.},
			KmerFreqMap{"A": 0.25, "G": 0.75},
			nil,
			errors.New("total freq for all kmers is not 1"),
		},
		// If the Input library contains a non-ACGT kmer
		{KmerFreqMap{"A": 0.5, "C": 0.5},
			KmerFreqMap{"A": 0.25, "N": 0.75},
			nil,
			errors.New("non-ACGT kmer"),
		},
		// If the Input library contains different length kmers
		{KmerFreqMap{"A": 0.5, "C": 0.5},
			KmerFreqMap{"A": 0.25, "C": 0.75, "GT": 0.1},
			nil,
			errors.New("KmerFreqMap contains non-length"),
		},
	}

	for _, table := range tables {
		actual, err := KmerEnrichments(table.pdMap, table.inputMap)
		// If the expected error is not nil,
		if table.error != nil {
			assert.Error(t, err)
			assert.Contains(t, err.Error(), table.error.Error())
		} else {
			if !reflect.DeepEqual(actual, table.rMap) {
				t.Errorf("Enrichments from PD = %v; Input = %v incorrect, got: %v, want: %v.",
					table.pdMap, table.inputMap, actual, table.rMap)
			}
		}
	}
}
