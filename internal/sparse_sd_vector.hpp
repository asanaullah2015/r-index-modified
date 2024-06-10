// Copyright (c) 2017, Nicola Prezza.  All rights reserved.
// Use of this source code is governed
// by a MIT license that can be found in the LICENSE file.

/*
  * sparse_sd_vector: a wrapper on sd_vector of the sdsl library, with support for rank/select1
  */

//============================================================================


#ifndef INTERNAL_SPARSE_SD_VECTOR_HPP_
#define INTERNAL_SPARSE_SD_VECTOR_HPP_

#include <vector>

//using namespace std;
//using namespace sdsl;

#ifndef ulint
typedef uint64_t ulint;
#endif

#ifndef uint
typedef uint32_t uint;
#endif

namespace ri{

class sparse_sd_vector{

public:

	/*
	 * empty constructor. Initialize bitvector with length 0.
	 */
	sparse_sd_vector(){}

	/*
	 * constructor. build bitvector given a vector of bools
	 */
	sparse_sd_vector(std::vector<bool> &b){

		if(b.size()==0) return;

		u = b.size();

		sdsl::bit_vector bv(b.size());

		for(uint64_t i=0;i<b.size();++i)
			bv[i] = b[i];

		sdv = sdsl::sd_vector<>(bv);
		rank1 = sdsl::sd_vector<>::rank_1_type(&sdv);
		select1 = sdsl::sd_vector<>::select_1_type(&sdv);

	}

	/*
	 * constructor. build bitvector given a bit_vector
	 */
	sparse_sd_vector(sdsl::bit_vector &bv){

		sdv = sdsl::sd_vector<>(bv);
		rank1 = sdsl::sd_vector<>::rank_1_type(&sdv);
		select1 = sdsl::sd_vector<>::select_1_type(&sdv);

	}

	sparse_sd_vector & operator= (const sparse_sd_vector & other) {

		u = other.sdv.size();
		sdv = sdsl::sd_vector<>(other.sdv);
		rank1 = sdsl::sd_vector<>::rank_1_type(&sdv);
		select1 = sdsl::sd_vector<>::select_1_type(&sdv);

	    return *this;
	}

	/*
	 * not implemented
	 * argument: a boolean b
	 * behavior: append b at the end of the bitvector.
	 */
	//void push_back(bool b){}

	/*
	 * argument: position i in the bitvector
	 * returns: bit in position i
	 * only access! the bitvector is static.
	 */
	bool operator[](const ulint i) const {

		assert(i<size());
		return sdv[i];

	}

	bool at(const ulint i) const {
		return operator[](i);
	}

	/*
	 * argument: position i in the bitvector
	 * returns: number of bits equal to 1 before position i excluded
	 */
	ulint rank(const ulint i) const {

		assert(i<=size());
		return rank1(i);

	}

	/*
	 * input: position 0<=i<=n
	 * output: predecessor of i (position i excluded) in
	 * rank space (apply select to get bitvector space)
	 */
	ulint predecessor_rank(const ulint i) const {

		/*
		 * i must have a predecessor
		 */
		assert(rank(i)>0);

		return rank(i)-1;

	}

	/*
	 * input: position 0<=i<=n
	 * output: predecessor of i (i excluded) in
	 * bitvector space
	 */
	ulint predecessor(const ulint i) const {

		/*
		 * i must have a predecessor
		 */
		assert(rank(i)>0);

		return select(rank(i)-1);


	}

	/*
	 * input: position 0<=i<=n
	 * output: rank of predecessor of i (i excluded) in
	 * bitvector space. If i does not have a predecessor,
	 * return rank of the last bit set in the bitvector
	 */
	ulint predecessor_rank_circular(const ulint i) const {

		return rank(i)==0 ? number_of_1()-1 : rank(i)-1;

	}

	/*
	 * retrieve length of the i-th gap (i>=0). gap length includes the leading 1
	 * \param i<number_of_1()
	 *
	 */
	ulint gapAt(const ulint i) const {

		assert(i<number_of_1());

		if(i==0) return select(0)+1;

		return select(i)-select(i-1);

	}

	/*
	 * argument: integer i>0
	 * returns: position of the i-th one in the bitvector. i starts from 0!
	 */
	ulint select(const ulint i) const {

		assert(i<number_of_1());
		return select1(i+1);//in sd_vector, i starts from 1

	}

	/*
	* returns: size of the bitvector
	*/
	ulint size() const {return u;}

	/*
	 * returns: number of 1s in the bitvector
	 */
	ulint number_of_1() const {return rank1(size()); }

	/* serialize the structure to the ostream
	 * \param out	 the ostream
	 */
	ulint serialize(std::ostream& out) const {

		ulint w_bytes = 0;

		out.write((char*)&u, sizeof(u));

		w_bytes += sizeof(u);

		if(u==0) return w_bytes;

		w_bytes += sdv.serialize(out);

		return w_bytes;

	}

	/* load the structure from the istream
	 * \param in the istream
	 */
	void load(std::istream& in) {

		in.read((char*)&u, sizeof(u));

		if(u==0) return;

		sdv.load(in);
		rank1 = sdsl::sd_vector<>::rank_1_type(&sdv);
		select1 = sdsl::sd_vector<>::select_1_type(&sdv);

	}

private:

	//bitvector length
	ulint u = 0;

	sdsl::sd_vector<> sdv;
	sdsl::sd_vector<>::rank_1_type rank1;
	sdsl::sd_vector<>::select_1_type select1;

};

}


#endif /* INTERNAL_SPARSE_SD_VECTOR_HPP_ */
