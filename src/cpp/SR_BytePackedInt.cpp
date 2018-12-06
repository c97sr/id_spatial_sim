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


#include"SR_BytePackedInt.h"
#include<limits>

SR::BytePackedInt::BytePackedInt() {
	baseInt = 0;
	intMax = numeric_limits<unsigned int>::max();
};

unsigned int SR::BytePackedInt::IntPow(unsigned int x, unsigned int a) {
	unsigned int rtn=1;
	for (unsigned int i=1;i<=a;i++) rtn *= x;
	return rtn;
};

unsigned int SR::BytePackedInt::Get(unsigned int l_bound, unsigned int u_bound){
#ifdef _DEBUG
	if (u_bound < l_bound || l_bound < 1 || u_bound > noBits) SR::srerror("Range problems in SR_BytePackedInt.");
#endif
	unsigned int rtn,mask;
	mask = IntPow(intTwo,u_bound)-IntPow(intTwo,l_bound-1);
	rtn = mask & baseInt;
	rtn /= IntPow(intTwo,l_bound-1);
	return rtn;
};

void SR::BytePackedInt::Set(unsigned int l_bound, unsigned int u_bound, unsigned int value){
#ifdef _DEBUG
	if (u_bound < l_bound || l_bound < 1) SR::srerror("Range problems in SR_BytePackedInt.");
	if (value > IntPow(intTwo,u_bound-l_bound+1) - 1) SR::srerror("Value too large in SR_BytePackedInt.");
#endif
	unsigned int mask;
	mask = IntPow(intTwo,u_bound)-IntPow(intTwo,l_bound-1);
	value = value*IntPow(intTwo,l_bound-1);
	baseInt = (intMax ^ mask) & (baseInt);
	baseInt += value;
};

void SR::BytePackedInt::Increment(unsigned int l_bound, unsigned int u_bound) {
	unsigned int tmp;
	tmp = Get(l_bound, u_bound);
	tmp++;
	Set(l_bound,u_bound,tmp);
};

void SR::BytePackedInt::Decrement(unsigned int l_bound, unsigned int u_bound) {
	unsigned int tmp;
	tmp = Get(l_bound, u_bound);
	tmp--;
	Set(l_bound,u_bound,tmp);
};

void SR::BytePackedInt::FastIncrement(unsigned int l_bound, unsigned int u_bound) {
	baseInt+=IntPow(intTwo,l_bound-1);
};

void SR::BytePackedInt::FastDecrement(unsigned int l_bound, unsigned int u_bound) {
	baseInt-=IntPow(intTwo,l_bound-1);
};

unsigned int SR::BytePackedInt::GetBasicInt() {
	return baseInt;
};

void SR::BytePackedInt::SetBasicInt(unsigned int bii) {
	baseInt=bii;
};
