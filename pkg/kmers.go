package rbns

import (
	"fmt"
	//"math"
)

type KmerCountMap map[string]int
type KmerFreqMap map[string]float64
type KmersSlice []string


func AllKmers(k int) KmersSlice {
	if k < 1 {
		panic("k for AllKmers must be >= 1")
	}
	//nKmers := int(math.Pow(4, float64(k)))

	currKmers := KmersSlice{""}
	for i := 0; i < k; i++ {
		var kPlus1mers KmersSlice
		for _, nt := range []string{"A", "C", "G", "T"} {
			for _, currKmer := range currKmers {
				kPlus1mers = append(kPlus1mers, nt + currKmer)
			}
		}
		currKmers = kPlus1mers
	}

	return currKmers
}


func CountKmers(s string, k int) KmerCountMap {
	kmerCounts := make(KmerCountMap)
	nKmers := len(s) - k + 1
	// If the string is shorter than k, there are no kmers
	if nKmers < 1 {
		return kmerCounts
	}

	for i := 0; i < nKmers; i++ {
		kmer := s[i:i+k]
		kmerCounts[kmer]++
	}
	return kmerCounts
}

func isACGT(s string) bool {
	acgt := []string{"A", "C", "G", "T"}
	for _, b := range acgt {
		if b == s {
			return true
		}
	}
	return false
}

func isACGTonlyKmer(s string) bool {
	// Return false if s is empty
	if len(s) == 0 {
		return false
	}
	for _, c := range s {
		if !isACGT(string(c)) {
			return false
		}
	}
	return true
}

func FilterToACGTkmers(m KmerCountMap) KmerCountMap {
	filtKmerCounts := make(KmerCountMap)
	for kmer, count := range m {
		if isACGTonlyKmer(kmer) {
			filtKmerCounts[kmer] = count
		}
	}
	return filtKmerCounts
}


func CountsMapToFreqMap(m KmerCountMap) (KmerFreqMap, error) {
	// First get the total number of counts in the map
	var totalCounts int
	for kmer, count := range m {
		if !isACGTonlyKmer(kmer) {
			return nil, fmt.Errorf("KmerCountMap contains non-ACGT kmer %s", kmer)
		}
		if count < 0 {
			return nil, fmt.Errorf("KmerCountMap contains negative counts for %s", kmer)
		}
		totalCounts += count
	}

	kmerFreqMap := make(KmerFreqMap)
	for kmer, count := range m {
		kmerFreqMap[kmer] = float64(count) / float64(totalCounts)
	}
	return kmerFreqMap, nil
}
