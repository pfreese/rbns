package rbns

import (
	"reflect"
	"testing"
)

func TestCountKmers(t *testing.T) {
	tables := []struct {
		s string
		k int
		kmerCounts map[string]int
	}{
		{"ABCDE", 2, map[string]int{"AB": 1, "BC": 1, "CD": 1, "DE": 1}},
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