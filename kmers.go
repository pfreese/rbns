package rbns

import (
	"fmt"
	"math"
)

type KmerCountMap map[string]int
type KmerFreqMap map[string]float64
type KmerRMap map[string]float64
type KmersSlice []string


func AllKmers(k int) KmersSlice {
	if k < 1 {
		panic("k for AllKmers must be >= 1")
	}

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

func getk(m KmerFreqMap) int {
	for kmer, _ := range m {
		return len(kmer)
	}
	return 0
}

// checkValidKmerFreqMap performs 2 checks:
// 1. All key entries in m are of length k
// 2. All value (frequency) entries in m are non-negative and
//    together sum to 1.
func checkValidKmerFreqMap(m KmerFreqMap, k int) error {
	var totalFreq float64
	for kmer, freq := range m {
		if len(kmer) != k {
			return fmt.Errorf("KmerFreqMap contains non-length %d kmer: %s", k, kmer)
		}
		if !isACGTonlyKmer(kmer) {
			return fmt.Errorf("non-ACGT kmer: %s", kmer)
		}
		if freq < 0 {
			return fmt.Errorf("freq for %s is <0: %f", kmer, freq)
		}
		totalFreq += freq
	}
	// Make sure the frequencies sum to 1
	if math.Abs(totalFreq - 1.) > 1e-6 {
		return fmt.Errorf("total freq for all kmers is not 1 (=%f)", totalFreq)
	}
	return nil
}

func KmerEnrichments(pd, input KmerFreqMap) (KmerRMap, error) {
	// Get the kmer size
	k := getk(pd)

	err := checkValidKmerFreqMap(pd, k)
	if err != nil {
		return nil, err
	}
	err = checkValidKmerFreqMap(input, k)
	if err != nil {
		return nil, err
	}

	kmerRMap := make(KmerRMap)

	// Go through all kmers, determining if each is present in both
	// maps, and if so, calculating R.
	for _, kmer := range AllKmers(k) {

		pdFreq, pdOK := pd[kmer]
		inputFreq, inputOK := input[kmer]

		// If the kmer's absent in either PD or Input,
		// or if it has 0 frequency in the input, assign it R = 1
		if (!pdOK || !inputOK) || inputFreq == 0. {
			kmerRMap[kmer] = 1
		} else {
			// Expected R definition for normal cases
			kmerRMap[kmer] = pdFreq / inputFreq
		}
	}
	return kmerRMap, nil
}
