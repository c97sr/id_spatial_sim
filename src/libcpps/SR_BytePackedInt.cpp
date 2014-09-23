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
