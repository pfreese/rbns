package rbns

func CountKmers(s string, k int) map[string]int {
	kmerCounts := make(map[string]int)
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

func FilterToACGTkmers(m map[string]int) map[string]int {
	filtKmerCounts := make(map[string]int)
	for kmer, count := range m {
		if isACGTonlyKmer(kmer) {
			filtKmerCounts[kmer] = count
		}
	}
	return filtKmerCounts
}
