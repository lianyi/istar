#include <iostream>
#include <openbabel/obconversion.h>
#include <openbabel/mol.h>
using namespace std;
using namespace OpenBabel;

int main(int argc, char* argv[])
{
	if (argc < 2)
	{
		cout << "subset ZINC00537755.mol2" << endl;
		return 0;
	}

	const auto SmartsPatterns =
	{
		"[!#1]", // heavy
		"[#6+0!$(*~[#7,#8,F]),SH0+0v2,s+0,S^3,Cl+0,Br+0,I+0]", // hydrophobic
		"[a]", // aromatic
		"[$([O,S;H1;v2]-[!$(*=[O,N,P,S])]),$([O,S;H0;v2]),$([O,S;-]),$([N&v3;H1,H2]-[!$(*=[O,N,P,S])]),$([N;v3;H0]),$([n,o,s;+0]),F]", // acceptor
		"[N!H0v3,N!H0+v4,OH+0,SH+0,nH+0]", // donor
	};

	OBConversion obConversion;
	obConversion.SetInFormat(argv[1] + string(argv[1]).find_last_of('.') + 1);
	OBMol obMol;
	obConversion.ReadFile(&obMol, argv[1]);
	for (const auto& p : SmartsPatterns)
	{
		OBSmartsPattern smarts;
		smarts.Init(p);
		smarts.Match(obMol);
		for (const auto& map : smarts.GetMapList())
		{
			cout << map.front() << endl;
		}
		cout << endl;
	}
}
