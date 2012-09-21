#include "smiley.h"
#include "molecule.h"

#include "test.h"
#include <cassert>

using namespace Smiley;

MoleculeSmilesCallback callback;
const Molecule &mol = callback.molecule;

void parse(const std::string &str)
{
  std::cout << "Parsing: " << str << std::endl;
  Parser<MoleculeSmilesCallback> parser(callback);

  try {
    parser.parse(str);
  } catch (Exception &e) {
    if (e.type() == Exception::SyntaxError)
      std::cout << "Syntax";
    else
      std::cout << "Semantics";
    std::cout << "Error: " << e.what() << "." << std::endl;
    assert(e.pos() < str.size());
    assert(e.length() < e.pos() + str.size());
    std::cout << str << std::endl;
    for (std::size_t i = 0; i < e.pos(); ++i)
      std::cout << " ";
    // print error indicater
    for (std::size_t i = 0; i < e.length(); ++i)
      std::cout << "^";
    std::cout << std::endl;
    // print positions for tens
    for (std::size_t i = 0; i < str.size(); ++i)
      if ((i % 10) == 0)
        std::cout << i / 10;
      else
        std::cout << " ";
    std::cout << std::endl;
    // print positions for units
    for (std::size_t i = 0; i < str.size(); ++i)
      std::cout << i % 10;
    std::cout << std::endl;
  }
}

void testException(const std::string &str, ErrorCode errorCode, bool expect = true)
{
  std::cout << "Parsing: " << str << std::endl;
  MoleculeSmilesCallback callback;
  Parser<MoleculeSmilesCallback> parser(callback);

  if (!expect)
    parser.disableExceptions(errorCode);

  try {
    parser.parse(str);
    if (expect)
      ASSERT(std::string("Expected Exception not thrown").empty());
  } catch (const Exception &e) {
    if (!expect)
      ASSERT(std::string("Unexpected Exception thrown").empty());
    ASSERT(e.errorCode() == errorCode);
  }
}

void parseElement(const std::string &str, int element, bool aromatic = false)
{
  parse(str);
  COMPARE(mol.atoms[0].element, element);
  COMPARE(mol.atoms[0].aromatic, aromatic);
}

void parseBond(const std::string &str, int order)
{
  parse(str);
  COMPARE(mol.bonds.size(), 1);
  COMPARE(mol.bonds[0].order, order);
}

int main(int argc, char **argv)
{
  //////////////////////////////////////////////
  //
  //
  //  Test Syntax
  //
  //
  //////////////////////////////////////////////

  parse("");

  parse("C");
  COMPARE(mol.atoms.size(), 1);
  COMPARE(mol.bonds.size(), 0);
  COMPARE(mol.atoms[0].element, 6);
  // check for default organic subset values once
  COMPARE(mol.atoms[0].isotope, -1);
  COMPARE(mol.atoms[0].aromatic, false);
  COMPARE(mol.atoms[0].hCount, -1);
  COMPARE(mol.atoms[0].charge, 0);
  COMPARE(mol.atoms[0].atomClass, 0);

  parse("C ");
  parse("C\t");
  parse("C\n");
  parse("C\r");

  parse("*");
  COMPARE(mol.atoms.size(), 1);
  COMPARE(mol.atoms[0].element, 0);

  parse("CC");
  COMPARE(mol.atoms.size(), 2);
  COMPARE(mol.bonds.size(), 1);
  ASSERT(mol.bond(0, 1) != 0);
  COMPARE(mol.bond(0, 1)->order, 1);

  // syntax error: no symbol
  parse("[]");

  parse("[C]");
  COMPARE(mol.atoms.size(), 1);
  COMPARE(mol.bonds.size(), 0);
  COMPARE(mol.atoms[0].element, 6);
  // check for default organic subset values once
  COMPARE(mol.atoms[0].isotope, -1);
  COMPARE(mol.atoms[0].aromatic, false);
  COMPARE(mol.atoms[0].hCount, 0);
  COMPARE(mol.atoms[0].charge, 0);
  COMPARE(mol.atoms[0].atomClass, 0);

  // syntax error: no symbol
  parse("[13]");
  COMPARE(mol.atoms.size(), 0);

  parse("[13C]");
  COMPARE(mol.atoms[0].element, 6);
  COMPARE(mol.atoms[0].isotope, 13);
  parse("[013C]");
  COMPARE(mol.atoms[0].element, 6);
  COMPARE(mol.atoms[0].isotope, 13);
  parse("[000C]");
  COMPARE(mol.atoms[0].element, 6);
  COMPARE(mol.atoms[0].isotope, 0);
  parse("[999C]");
  COMPARE(mol.atoms[0].element, 6);
  COMPARE(mol.atoms[0].isotope, 999);

  parse("[C@](*)(*)(*)*");
  COMPARE(mol.chirality.begin()->second.first, AntiClockwise);
  parse("[C@@](*)(*)(*)*");
  COMPARE(mol.chirality.begin()->second.first, Clockwise);
  parse("[C@TH1](*)(*)(*)*");
  COMPARE(mol.chirality.begin()->second.first, TH1);
  parse("[C@TH2](*)(*)(*)*");
  COMPARE(mol.chirality.begin()->second.first, TH2);
  parse("FC(Cl)=[C@AL1]=C(Br)I");
  COMPARE(mol.chirality.begin()->second.first, AL1);
  parse("FC(Cl)=[C@AL2]=C(Br)I");
  COMPARE(mol.chirality.begin()->second.first, AL2);
  parse("[C@SP1](*)(*)(*)*");
  COMPARE(mol.chirality.begin()->second.first, SP1);
  parse("[C@SP2](*)(*)(*)*");
  COMPARE(mol.chirality.begin()->second.first, SP2);
  parse("[C@SP3](*)(*)(*)*");
  COMPARE(mol.chirality.begin()->second.first, SP3);
  parse("[C@TB1](*)(*)(*)(*)*");
  COMPARE(mol.chirality.begin()->second.first, TB1);
  parse("[C@TB10](*)(*)(*)(*)*");
  COMPARE(mol.chirality.begin()->second.first, TB10);
  parse("[C@TB20](*)(*)(*)(*)*");
  COMPARE(mol.chirality.begin()->second.first, TB20);
  parse("[C@OH1](*)(*)(*)(*)(*)*");
  COMPARE(mol.chirality.begin()->second.first, OH1);
  parse("[C@OH10](*)(*)(*)(*)(*)*");
  COMPARE(mol.chirality.begin()->second.first, OH10);
  parse("[C@OH20](*)(*)(*)(*)(*)*");
  COMPARE(mol.chirality.begin()->second.first, OH20);
  parse("[C@OH30](*)(*)(*)(*)(*)*");
  COMPARE(mol.chirality.begin()->second.first, OH30);

  parse("[C]");
  COMPARE(mol.atoms[0].hCount, 0);
  parse("[CH]");
  COMPARE(mol.atoms[0].hCount, 1);
  parse("[CH2]");
  COMPARE(mol.atoms[0].hCount, 2);
  parse("[CH4]");
  COMPARE(mol.atoms[0].hCount, 4);

  parse("[C-]");
  COMPARE(mol.atoms[0].charge, -1);
  parse("[C--]");
  COMPARE(mol.atoms[0].charge, -2);
  parse("[C-1]");
  COMPARE(mol.atoms[0].charge, -1);
  parse("[C-3]");
  COMPARE(mol.atoms[0].charge, -3);
  parse("[C+]");
  COMPARE(mol.atoms[0].charge, 1);
  parse("[C++]");
  COMPARE(mol.atoms[0].charge, 2);
  parse("[C+1]");
  COMPARE(mol.atoms[0].charge, 1);
  parse("[C+2]");
  COMPARE(mol.atoms[0].charge, 2);

  parse("[C:1]");
  COMPARE(mol.atoms[0].atomClass, 1);
  parse("[C:001]");
  COMPARE(mol.atoms[0].atomClass, 1);
  parse("[C:123]");
  COMPARE(mol.atoms[0].atomClass, 123);
  parse("[C:999]");
  COMPARE(mol.atoms[0].atomClass, 999);

  parseElement("B", B);
  parseElement("C", C);
  parseElement("N", N);
  parseElement("O", O);
  parseElement("P", P);
  parseElement("S", S);
  parseElement("F", F);
  parseElement("Cl", Cl);
  parseElement("Br", Br);
  parseElement("I", I);
  parseElement("b", B, true);
  parseElement("c", C, true);
  parseElement("n", N, true);
  parseElement("o", O, true);
  parseElement("p", P, true);
  parseElement("s", S, true);
  parseElement("[H]", H);
  parseElement("[Li]", Li);
  parseElement("[Na]", Na);
  parseElement("[Xe]", Xe);
  parseElement("[Cu]", Cu);
  parseElement("[Y]", Y);
  parseElement("[W]", W);
  parseElement("[Cd]", Cd);
  parseElement("[Rh]", Rh);
  parseElement("[Ga]", Ga);
  parseElement("[Ba]", Ba);
  parseElement("[F]", F);
  parseElement("[S]", S);
  parseElement("[c]", C, true);
  parseElement("[n]", N, true);
  parseElement("[o]", O, true);
  parseElement("[p]", P, true);
  parseElement("[s]", S, true);
  parseElement("[se]", Se, true);
  parseElement("[as]", As, true);
  parseElement("[*]", 0);

  parseBond("**", 1);
  parseBond("*-*", 1);
  parseBond("*=*", 2);
  parseBond("*#*", 3);
  parseBond("*$*", 4);
  parseBond("*:*", 5);


  parse("C1CC1");
  COMPARE(mol.bonds.size(), 3);
  COMPARE(mol.bonds[0].source, 0);
  COMPARE(mol.bonds[0].target, 1);
  COMPARE(mol.bonds[0].order, 1);
  COMPARE(mol.bonds[1].source, 1);
  COMPARE(mol.bonds[1].target, 2);
  COMPARE(mol.bonds[1].order, 1);
  COMPARE(mol.bonds[2].source, 0);
  COMPARE(mol.bonds[2].target, 2);
  COMPARE(mol.bonds[2].order, 1);

  parse("C=1CC1");
  COMPARE(mol.bonds.size(), 3);
  COMPARE(mol.bonds[2].order, 2);
  
  parse("C1CC=1");
  COMPARE(mol.bonds.size(), 3);
  COMPARE(mol.bonds[2].order, 2);

  parse("C1C2CCCC2C1");
  COMPARE(mol.bonds.size(), 8);
  ASSERT(mol.bond(0, 6) != 0);
  ASSERT(mol.bond(1, 5) != 0);

  parse("CC(C)C");
  COMPARE(mol.atoms.size(), 4);
  COMPARE(mol.bonds.size(), 3);
  ASSERT(mol.bond(0, 1) != 0);
  ASSERT(mol.bond(1, 2) != 0);
  ASSERT(mol.bond(1, 3) != 0);

  parse("CC(CC");
  parse("CC)CC");
  parse("C1CCC");
  parse("C-1CCCCC=1");

  parse("CC(C)(C)C");
  COMPARE(mol.atoms.size(), 5);
  COMPARE(mol.bonds.size(), 4);
  ASSERT(mol.bond(0, 1) != 0);
  ASSERT(mol.bond(1, 2) != 0);
  ASSERT(mol.bond(1, 3) != 0);
  ASSERT(mol.bond(1, 4) != 0);

  parse("C(C(C(C)))C");
  COMPARE(mol.atoms.size(), 5);
  COMPARE(mol.bonds.size(), 4);
  ASSERT(mol.bond(0, 1) != 0);
  ASSERT(mol.bond(1, 2) != 0);
  ASSERT(mol.bond(2, 3) != 0);
  ASSERT(mol.bond(0, 4) != 0);

  parse("C[C]C");
  COMPARE(mol.atoms.size(), 3);
  parse("C[C](C)C");
  COMPARE(mol.atoms.size(), 4);

  
  //parse("[Cl-].COC(=O)C12[NH2+]CC3(C)C(CCC23C)C1 10020");
  //return 0;

  //////////////////////////////////////////////
  //
  //
  //  Test Semantics
  //
  //
  //////////////////////////////////////////////

  //
  // Tetrahedral
  //
  std::map<int, std::pair<Chirality, std::vector<int> > >::const_iterator chiralInfo;
  parse("C[C@](F)(Cl)Br");
  COMPARE(mol.chirality.size(), 1);
  chiralInfo = mol.chirality.find(1);
  ASSERT(chiralInfo != mol.chirality.end());
  COMPARE(chiralInfo->second.first, AntiClockwise);
  COMPARE(chiralInfo->second.second.size(), 4);
  COMPARE(chiralInfo->second.second[0], 0);
  COMPARE(chiralInfo->second.second[1], 2);
  COMPARE(chiralInfo->second.second[2], 3);
  COMPARE(chiralInfo->second.second[3], 4);

  parse("C[C@H](F)Cl");
  COMPARE(mol.chirality.size(), 1);
  chiralInfo = mol.chirality.find(1);
  ASSERT(chiralInfo != mol.chirality.end());
  COMPARE(chiralInfo->second.first, AntiClockwise);
  COMPARE(chiralInfo->second.second.size(), 4);
  COMPARE(chiralInfo->second.second[0], 0);
  COMPARE(chiralInfo->second.second[1], implicitHydrogen());
  COMPARE(chiralInfo->second.second[2], 2);
  COMPARE(chiralInfo->second.second[3], 3);

  parse("C[C@TH1](F)(Cl)Br");
  COMPARE(mol.chirality.size(), 1);
  chiralInfo = mol.chirality.find(1);
  ASSERT(chiralInfo != mol.chirality.end());
  COMPARE(chiralInfo->second.first, TH1);
  COMPARE(chiralInfo->second.second.size(), 4);
 
  parse("C[C@TH2](F)(Cl)Br");
  COMPARE(mol.chirality.size(), 1);
  chiralInfo = mol.chirality.find(1);
  ASSERT(chiralInfo != mol.chirality.end());
  COMPARE(chiralInfo->second.first, TH2);
  COMPARE(chiralInfo->second.second.size(), 4);

  //
  // Allene
  //
  parse("NC(Cl)=[C@AL1]=C(F)Br");
  COMPARE(mol.chirality.size(), 1);
  chiralInfo = mol.chirality.find(3);
  ASSERT(chiralInfo != mol.chirality.end());
  COMPARE(chiralInfo->second.first, AL1);
  COMPARE(chiralInfo->second.second.size(), 4);
  COMPARE(chiralInfo->second.second[0], 0);
  COMPARE(chiralInfo->second.second[1], 2);
  COMPARE(chiralInfo->second.second[2], 5);
  COMPARE(chiralInfo->second.second[3], 6);
  
  parse("NC(Cl)=[C@]=C(F)Br");
  parse("NC(Cl)=[C@@]=C(F)Br");

  //
  // Square Planar
  //
  parse("C[Pt@SP1](F)(Cl)Br");
  COMPARE(mol.chirality.size(), 1);
  chiralInfo = mol.chirality.find(1);
  ASSERT(chiralInfo != mol.chirality.end());
  COMPARE(chiralInfo->second.first, SP1);
  COMPARE(chiralInfo->second.second.size(), 4);

  //
  // Trigonal Bypiramidal
  //
  parse("*[*@TB7](*)(*)(*)*");
  COMPARE(mol.chirality.size(), 1);
  chiralInfo = mol.chirality.find(1);
  ASSERT(chiralInfo != mol.chirality.end());
  COMPARE(chiralInfo->second.first, TB7);
  COMPARE(chiralInfo->second.second.size(), 5);

  //
  // Octahedral
  //
  parse("*[*@OH23](*)(*)(*)(*)*");
  COMPARE(mol.chirality.size(), 1);
  chiralInfo = mol.chirality.find(1);
  ASSERT(chiralInfo != mol.chirality.end());
  COMPARE(chiralInfo->second.first, OH23);
  COMPARE(chiralInfo->second.second.size(), 6);


  parse("C[C@H]1(F)CCCC1");
  parse("C1CCC[C@H]1F");
  parse("C[C@H]1(F)CCCCC[C@H]1Cl");

  parse("O=C(Cc1ccccc1F)N1CCC[C@H]1Cn1nc2CCCCCc2cc1=O 315718");
  parse("O=C(C1CCCN1)N1C[C@H]2C(=O)Nc3ccccc3C(=O)N[C@@]2(C)C1 356474");

  //
  // Test exceptions
  //
  testException("[HH1]", HydrogenHydrogenCount);
  testException("C1CC", UnmatchedRingBond);
  testException("C-1CCCCC=1", ConflictingRingBonds);
  testException("C12CCCCC12", InvalidRingBond);
  testException("C11", InvalidRingBond);
  testException("C[C@](F)Cl", InvalidChiralValence);
  testException("C[C@H2](F)Cl", InvalidChiralHydrogenCount);

  testException("C1CC", UnmatchedRingBond, false);
  testException("C-1CCCCC=1", ConflictingRingBonds, false);
  testException("C12CCCCC12", InvalidRingBond, false);
  testException("C11", InvalidRingBond, false);
  testException("C[C@](F)Cl", InvalidChiralValence, false);
  testException("C[C@H2](F)Cl", InvalidChiralHydrogenCount, false);

  parse("[F-][Ti+4]123456789([F-])C%10=C1[C]15=[C]8(CCCC1)[C-]4%10CC[C-]16C2=C3[C]27=[C]91CCCC2 36749863");
}
