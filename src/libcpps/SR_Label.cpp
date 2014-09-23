#include"SR_Label.h"

void SR::Label::Set(string s) {
	size = static_cast<int>(s.size());
	if (size > noChars) SR::srerror("String assigned to label is too long in SR::Label::Set(string s)");
	for (int i=0;i<size;++i)
		data[i]=s[i];
}

string SR::Label::Get() {
	string s(data.begin(),data.begin()+size);
	return s;
}

bool SR::Label::IsEqual(string s) {
	if (static_cast<int>(s.size()) != size) return false;
	for (int i=0;i<size;++i) if (s[i] != data[i]) return false;
	return true;
};

ofstream& SR::operator<<(ofstream& ofs, SR::Label& lab) {
	static char *filePointer;
	filePointer = (char*)(&lab.size); ofs.write(filePointer,sizeof(int));
	filePointer = (char*)(&lab.data[0]); ofs.write(filePointer,lab.noChars*sizeof(char));
	return ofs;
};

ifstream& SR::operator>>(ifstream& ifs, SR::Label& lab) {
	static char *filePointer;
	filePointer = (char*)(&lab.size); ifs.read(filePointer,sizeof(int));
	filePointer = (char*)(&lab.data[0]); ifs.read(filePointer,lab.noChars*sizeof(char));
	return ifs;
};
