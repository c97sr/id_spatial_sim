/*  Copyright 2018 Steven Riley.

    This file is part of id_spatial_sim.

    id_spatial_sim is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    id_spatial_sim is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with id_spatial_sim.  If not, see <https://www.gnu.org/licenses/>.  */


#ifndef SR_INC_BYTEPACKEDINT
#define SR_INC_BYTEPACKEDINT

#include<iostream>
#include<limits>
// #include"nr.h"
#include"SR_Utility.h"

/*
 This is a self contained class within which ints can be aggregated to save memory.
 Note that each call needs to use the correct l_bound and u_bound - no record of where
 the numbers "are" is kept.  Any such record would add to the size of the class and thus
 negate a lot of the reason for having it!

 Usually will be used with compiler option /D SR_BYTEPACKED so as to leave a faster, larger
 version of the code.

 S Riley 26 July 2004

 */


using namespace std;

namespace SR {
	class BytePackedInt {
	private:
		unsigned int baseInt;
		unsigned int intMax;  				// Now set in the constructor
		static const unsigned int intTwo = 2;
		static const unsigned int noBits = 32;			// Change this if move to 64 bit machine
		unsigned int IntPow(unsigned int x, unsigned int a);
	public:
		BytePackedInt();
		unsigned int Get(unsigned int l_bound, unsigned int u_bound);
		void Set(unsigned int l_bound, unsigned int u_bound, unsigned int value);
		void Increment(unsigned int l_bound, unsigned int u_bound);
		void Decrement(unsigned int l_bound, unsigned int u_bound);
		void FastIncrement(unsigned int l_bound, unsigned int u_bound);
		void FastDecrement(unsigned int l_bound, unsigned int u_bound);
		unsigned int GetBasicInt();
		void SetBasicInt(unsigned int bii);
	};
}

#endif

