#include<iostream>

using namespace std;

typedef class B;

class A {
	friend class B;
	private:
		static B* x;
};

class B {
	friend class A;
	private:
		static A* y;
};

int main() {
	cout << "Hello world!" << endl;
	return 0;
}
