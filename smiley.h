/**
 * Copyright (c) 2012, Tim Vandermeersch
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the <organization> nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
#ifndef SMILEY_SMILEY_H
#define SMILEY_SMILEY_H

#include <string>
#include <vector>
#include <map>
#include <iostream>

//#include <cassert>
#define DEBUG 0

namespace Smiley {

  // OpenSMILES comments:
  // * chiral: TB only has 20 values (not 30) [OpenSMILES 2.2]
  // * @ is anti-clockwise and @@ is clockwise [OpenSMILES 3.9.6]
  //   exchange 'clockwise' and 'anti-clockwise' in 2nd sentence

  /**
   * @mainpage 
   *
   * Smiley is a SMILES/SMARTS parser that is OpenSMILES (http://www.opensmiles.org)
   * compliant (draft 2007-11-13).
   *
   * @section usage Usage
   *
   * Smiley is a single header file which contains the Parser class in the
   * Smiley namespace. User actions are invoked using a functor object that is
   * specified when a Parser is instantiated. Below is a simple example which
   * uses the PrintCallback to print the events as they occur. The example below
   * also shows how to handle exceptions and print useful error reports.
   *
   * @code
   * #include "smiley.h"
   * 
   * using namespace Smiley;
   * 
   * int main(int argc, char **argv)
   * {
   *   PrintCallback callback;
   *   Parser<PrintCallback> parser(callback);
   * 
   *   try {
   *     parser.parse(argv[1]);
   *   } catch (Exception &e) {
   *     if (e.type() == Exception::SyntaxError)
   *       std::cerr << "Syntax";
   *     else
   *       std::cerr << "Semantics";
   *     std::cerr << "Error: " << e.what() << "." << std::endl;
   *     std::cerr << argv[1] << std::endl;
   *     for (std::size_t i = 0; i < e.pos(); ++i)
   *       std::cerr << " ";
   *     for (std::size_t i = 0; i < e.length(); ++i)
   *       std::cerr << "^";
   *     std::cerr << std::endl;
   *   }
   * }
   * @endcode
   *
   * Example of errors:
   * @code
   * SyntaxError: Bracket atom expression contains invalid trailing characters.
   * F.FB(F)F.[NH2+251][C@@H](CP(c1ccccc1)c1ccccc1)C(C)(C)C 31586112
   *                ^^
   * SyntaxError: Unmatched branch opening.
   * CC(CC
   *   ^^^
   * SyntaxError: Unmatched branch closing.
   * CC)CC
   * ^^^
   * SemanticsError: Unmatched ring bond.
   * C1CCC
   *  ^
   * SemanticsError: Conflicing ring bonds.
   * C-1CCCCC=1
   *          ^
   * @endcode
   *
   * Although this is useful for debugging, using Smiley in an application
   * requires a custom callback function object.
   *
   * @subsection callback Callback
   * The callback function object type is a template parameter to the Parser
   * class. The should be derived publicly from CallbackBase so unneeded methods
   * do not have to be specified. For SMILES, the class may contain 4 methods:
   *
   * @code
   * struct MyCustomCallback
   * {
   *   void clear()
   *   {
   *     // prepare for new SMILES
   *   }
   *
   *   void addAtom(int element, bool aromatic, int isotope, int hCount, int charge, int atomClass)
   *   {
   *     // invoked when an atom is completely parsed
   *
   *     // element: 0, 1, 2, ... (0 is '*' in SMILES)
   *     // aromatic: true if the atom was lower case c, n, ...
   *     // isotope: -1, 0, 1, 2, 3, 4 ... (-1 means no isotope specified and is not the same as 0)
   *     // hCount: -1, 0, 1, 2, ..., 9 (-1 means default valence and is only for unbracketed atoms)
   *     // charge: -9, -8, ..., -1, 0, 1, ..., 8, 9
   *     // atomClass: 0, 1, 2, 3, 4, ... (0 means no atom class, specified atom classes should start from 1)
   *   }
   *
   *   void addBond(int source, int target, int order, bool isUp, bool isDown)
   *   {
   *     // invoked for each bond once both of it's atoms have been added by
   *     // calling addAtom(). This ensures that the bond's atom indexes are always valid.
   *
   *     // source: source atom index starting from 0 (order from addAtom() calls)
   *     // target: target atom index starting from 0 (order from addAtom() calls)
   *     // order: 1, 2, 3, 4, 5 (5 means aromatic)
   *     // isUp: true if bond is single order up bond '/'
   *     // isDown: true if bond is single order down bond '\'
   *   }
   *
   *   void setChiral(int index, Chirality chirality, const std::vector<int> &chiralNbrs)
   *   {
   *     // invoked at the end of parsing for each chiral atom
   *
   *     // index: atom index starting from 0
   *     // chirality: Clockwise, AntiClockwise, TH1, AL2, SP3, TB14, OH26, ...
   *     // chiralNbrs: atom indices of neighbors, size 4 for TH, AL and SP, size 5 for TB and 6 for OH
   *   }
   * };
   * @endcode
   *
   * @section smiles_grammar SMILES Syntax and Grammar
   *
   * Notations:
   * @code
   * TYPES
   *
   *   identifier   name for definition
   *   '...'        literal string
   *   expr         expression (may contain operators)
   *
   * STATEMENTS
   *
   *   identifier ::= definition (expressions)
   *
   * OPERATORS
   *
   *   expr expr      concatination is implicit
   *   expr | expr    or
   *   expr?          none or one
   *   expr+          one or more
   *   expr*          none or more
   * @endcode
   *
   * Grammar:
   * @code
   * PRIMITIVES
   *
   *   DIGIT ::= '0' | '1' | '2' | '3' | '4' | '5' | '6' | '7' | '8' | '9'
   *
   *   NUMBER ::= DIGIT+
   *
   * ATOMS
   *
   *   atom ::= bracket_atom | aliphatic_organic | aromatic_organic | '*'
   *
   * ORGANIC SUBSET ATOMS
   *
   *   aliphatic_organic ::= 'B' | 'C' | 'N' | 'O' | 'S' | 'P' | 'F' | 'Cl' | 'Br' | 'I'
   *   aromatic_organic ::= 'b' | 'c' | 'n' | 'o' | 's' | 'p'
   *
   * BRACKET ATOMS
   *
   *   bracket_atom ::= '[' isotope? symbol chiral? hcount? charge? class? ']'
   *   symbol := element_symbols | aromatic_symbols | '*'
   *   isotope ::= NUMBER
   *   element_symbols ::= 'H' | 'He' | 'Li' | 'Be' | 'B' | 'C' | 'N' | 'O' | 'F'
   *       | 'Ne' | 'Na' | 'Mg' | 'Al' | 'Si' | 'P' | 'S' | 'Cl' | 'Ar' | 'K'
   *       | 'Ca' | 'Sc' | 'Ti' | 'V' | 'Cr' | 'Mn' | 'Fe' | 'Co' | 'Ni' | 'Cu'
   *       | 'Zn' | 'Ga' | 'Ge' | 'As' | 'Se' | 'Br' | 'Kr' | 'Rb' | 'Sr' | 'Y'
   *       | 'Zr' | 'Nb' | 'Mo' | 'Tc' | 'Ru' | 'Rh' | 'Pd' | 'Ag' | 'Cd' | 'In'
   *       | 'Sn' | 'Sb' | 'Te' | 'I' | 'Xe' | 'Cs' | 'Ba' | 'Hf' | 'Ta' | 'W'
   *       | 'Re' | 'Os' | 'Ir' | 'Pt' | 'Au' | 'Hg' | 'Tl' | 'Pb' | 'Bi' | 'Po'
   *       | 'At' | 'Rn' | 'Fr' | 'Ra' | 'Rf' | 'Db' | 'Sg' | 'Bh' | 'Hs' | 'Mt'
   *       | 'Ds' | 'Rg' | 'La' | 'Ce' | 'Pr' | 'Nd' | 'Pm' | 'Sm' | 'Eu' | 'Gd'
   *       | 'Tb' | 'Dy' | 'Ho' | 'Er' | 'Tm' | 'Yb' | 'Lu' | 'Ac' | 'Th' | 'Pa'
   *       | 'U' | 'Np' | 'Pu' | 'Am' | 'Cm' | 'Bk' | 'Cf' | 'Es' | 'Fm' | 'Md'
   *       | 'No' | 'Lr'
   *   aromatic_symbols ::= 'c' | 'n' | 'o' | 'p' | 's' | 'se' | 'as'
   *
   * CHIRALITY
   *
   *   chiral ::= '@' | '@@'
   *       | '@TH1' | '@TH2'
   *       | '@AL1' | '@AL2'
   *       | '@SP1' | '@SP2' | '@SP3'
   *       | '@TB1' | '@TB2' | '@TB3' | ... | '@TB19' | '@TB20'
   *       | '@OH1' | '@OH2' | '@OH3' | ... | '@OH29' | '@OH30'
   *
   * HYDROGENS
   *
   *   hcount ::= 'H' | 'H' DIGIT
   *
   * CHARGE
   *
   *   charge ::= '-' | '-' DIGIT | '+' | '+' DIGIT
   *       | '--' *deprecated*
   *       | '++' *deprecated*
   *
   * ATOM CLASS
   *
   *   class ::= ':' NUMBER
   *
   * BONDS AND CHAINS
   *
   *   bond ::= '-' | '=' | '#' | '$' | ':' | '/' | '\'
   *   ringbond ::= bond? DIGIT | bond? '%' DIGIT DIGIT
   *   branched_atom ::= atom ringbond* branch*
   *   branch ::= '(' chain ')' | '(' bond chain ')' | '(' dot chain ')'
   *   chain ::= branched_atom | chain branched_atom | chain bond branched_atom | chain dot branched_atom
   *   dot ::= '.'
   *
   * SMILES STRINGS
   *
   *   smiles ::= chain terminator
   *   terminator ::= SPACE TAB | LINEFEED | CARRIAGE_RETURN | END_OF_STRING
   * @endcode
   * 
   * @section smiles_semantics SMILES Semantics
   *
   * For a detailed description of the OpenSMILES semantics, the specification
   * should be consulted. Apart for syntactical and grammatical correctness,
   * Smiley als verifies some basic semantics. When a semantics error is found,
   * an Exception is thrown with the corresponding ErrorCode. Some of these
   * exceptions (UnmatchedRingBond, ConflictingRingBond, InvalidRingBond,
   * InvalidChiralValence and InvalidChiralHydrogenCount) can be disabled using
   * Parser::disableExceptions(). By default, all exceptions are enabled to
   * ensure OpenSMILES compliance. The effect of disabling these exceptions is
   * specified in subsequent subsections.
   *
   * @subsection semantics_HH Hydrogen with Hydrogen Count
   * Hydrogen atoms can not have a hydrogen count. Hydrogen bound to a hydrogen
   * atom should be specified by two bracket atom expressions. If such an
   * expression is found during parsing, an Exception with ErrorCode
   * HydrogenHydrogenCount is thrown.
   *
   * Eamples:
   * @code
   * [HH]        invalid
   * [HH1]       invalid (same as [HH]
   * [HH3]       invalid
   * [HH0]       valid (same as [H])
   * [H][H]      valid
   * @endcode
   *
   * @subsection semantics_unmatched_ringbond Unmatched Ring Bond
   * When there is an unmatched ring bond, an Exception with ErrorCode
   * UnmatchedRingBond will be thrown. Disabling this exception will
   * ignore the ring bond(s).
   *
   * Example:
   * @code
   * C1CCC
   * @endcode
   *
   * @subsection semantics_conflicing_ringbond Conflicting Ring Bonds
   * When the bond type for ring bonds are explicitly specified at both ends,
   * these should be the same. If not, an Exception with ErrorCode
   * ConflictingRingBonds is thrown. Disabling this exception will use the
   * first bond specification (i.e. a single bond in the example below).
   *
   * Example:
   * @code
   * C-1CCCCCC=1
   * @endcode
   *
   * @subsection semantics_invalid_ringbond Invalid Ring Bonds
   * There are two types of invalid ring bonds. The first is when two atoms both
   * have the same two ring bonds. This would mean adding a parallel edge in the
   * graph which is not allowed. The second type is similar but results in a
   * self-loop by having a ring bond number twice. An Exception with ErrorCode
   * InvalidRingBond is thrown when such a bond is encountered. Disabling this
   * exception will ignore the invalid ring bonds.
   *
   * Eamples:
   * @code
   * C12CCCC12      parallel bond
   * C11            self-loop bond
   * @endcode
   *
   * @subsection semantics_chiral_valence Invalid Chiral Valence
   * When an atom is specified as being chiral, it should have the correct
   * number of neighboring atoms (possibly including an implicit H inside the
   * bracket. An Exception with ErrorCode InvalidChiralValence is thrown when
   * a valence is incorrect. Disabling this exception will result in no
   * verification and leaves this to the Callback functor.
   *
   * The valid valences are:
   * @code
   * Tetrahedral (TH)          : 4
   * Allene (AL)               : 4 (*)
   * Square Planar (SP)        : 4
   * Trigonal Bypiramidal (TB) : 5
   * Octahedral(OH)            : 6
   *
   * (*) The chiral atom has only 2 bonds but the neighbor's neighbors are
   *     counted: NC(Br)=[C@AL1]=C(F)I
   * @endcode
   *
   * @subsection semantics_chiral_hydrogens Invalid Chiral Hydrogen Count
   * Chiral atoms can only have one hydrogen in their bracket since multiple
   * hydrogens would make them not chiral. An Exception with ErrorCode
   * InvalidChiralHydrogenCount is thrown when such an atom is encountered.
   * Disabling this exception allows a chiral specification with a hydrogen
   * count higher than 2. Note: only 1 implcitHydrogen() is added to nbrs.
   *
   * Example:
   * @code
   * C[C@H2]F
   * @endcode
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *  
   * @section smarts_grammar SMARTS Grammar
   *
   * @code
   * PRIMITIVES
   *
   *   DIGIT ::= '0' | '1' | '2' | '3' | '4' | '5' | '6' | '7' | '8' | '9'
   *
   *   NUMBER ::= DIGIT+
   *
   * ATOMS
   *
   *   atom ::= bracket_atom | aliphatic_organic | aromatic_organic | '*'
   *
   * ORGANIC SUBSET ATOMS
   *
   *   aliphatic_organic ::= 'B' | 'C' | 'N' | 'O' | 'S' | 'P' | 'F' | 'Cl' | 'Br' | 'I'
   *
   *   aromatic_organic ::= 'b' | 'c' | 'n' | 'o' | 's' | 'p'
   *
   * BRACKET ATOMS
   *
   *   bracket_atom ::= '[' atom_expression+ ']'
   *
   *   atom_expression ::= atom_primitive
   *       | recursive_smarts
   *       | unary_operator atom_primitive
   *       | atom_expression atom_primitive
   *       | atom_expression binary_operator atom_primitive
   *
   *   recursive_smarts ::= '$(' chain ')'
   *
   *   unary_operator ::= '!'
   *
   *   binary_operator ::= '&' | ';' | ','
   *
   * ATOM PRIMITIVES
   *
   *   atom_primitive ::= isotope | symbol | atomic_number | 'a' | 'A'
   *       | degree | valence | connectivity
   *       | total_hcount | implicit_hcount
   *       | ring_membership | ring_size | ring_connectivity 
   *       | charge | chiral | class
   *
   *   symbol := element_symbols | aromatic_symbols | '*'
   *      
   *   element_symbols ::= 'H' | 'He' | 'Li' | 'Be' | 'B' | 'C' | 'N' | 'O' | 'F'
   *       | 'Ne' | 'Na' | 'Mg' | 'Al' | 'Si' | 'P' | 'S' | 'Cl' | 'Ar' | 'K'
   *       | 'Ca' | 'Sc' | 'Ti' | 'V' | 'Cr' | 'Mn' | 'Fe' | 'Co' | 'Ni' | 'Cu'
   *       | 'Zn' | 'Ga' | 'Ge' | 'As' | 'Se' | 'Br' | 'Kr' | 'Rb' | 'Sr' | 'Y'
   *       | 'Zr' | 'Nb' | 'Mo' | 'Tc' | 'Ru' | 'Rh' | 'Pd' | 'Ag' | 'Cd' | 'In'
   *       | 'Sn' | 'Sb' | 'Te' | 'I' | 'Xe' | 'Cs' | 'Ba' | 'Hf' | 'Ta' | 'W'
   *       | 'Re' | 'Os' | 'Ir' | 'Pt' | 'Au' | 'Hg' | 'Tl' | 'Pb' | 'Bi' | 'Po'
   *       | 'At' | 'Rn' | 'Fr' | 'Ra' | 'Rf' | 'Db' | 'Sg' | 'Bh' | 'Hs' | 'Mt'
   *       | 'Ds' | 'Rg' | 'La' | 'Ce' | 'Pr' | 'Nd' | 'Pm' | 'Sm' | 'Eu' | 'Gd'
   *       | 'Tb' | 'Dy' | 'Ho' | 'Er' | 'Tm' | 'Yb' | 'Lu' | 'Ac' | 'Th' | 'Pa'
   *       | 'U' | 'Np' | 'Pu' | 'Am' | 'Cm' | 'Bk' | 'Cf' | 'Es' | 'Fm' | 'Md'
   *       | 'No' | 'Lr'
   *   
   *   aromatic_symbols ::= 'c' | 'n' | 'o' | 'p' | 's' | 'se' | 'as'
   *
   *   isotope ::= NUMBER
   *
   *   atomic_number ::= '#' NUMBER
   *
   *   degree ::= 'D' | 'D' DIGIT
   *
   *   valence ::= 'v' | 'v' NUMBER
   *
   *   connectivity ::= 'X' | 'X' NUMBER
   *
   *   total_hcount ::= 'H' | 'H' DIGIT
   *
   *   implicit_hcount ::= 'h' | 'h' DIGIT
   *
   *   ring_membership ::= 'R' | 'R' NUMBER
   *
   *   ring_size ::= 'r' | 'r' NUMBER
   *
   *   ring_connectivity ::= 'x' | 'x' NUMBER
   *
   *   charge ::= '-' | '-' DIGIT | '+' | '+' DIGIT
   *       | '--' *deprecated*
   *       | '++' *deprecated*
   *
   *   chiral ::= '@' | '@@'
   *       | '@TH1' | '@TH2'
   *       | '@AL1' | '@AL2'
   *       | '@SP1' | '@SP2' | '@SP3'
   *       | '@TB1' | '@TB2' | '@TB3' | ... | '@TB19' | '@TB20'
   *       | '@OH1' | '@OH2' | '@OH3' | ... | '@OH29' | '@OH30'
   *
   *   class ::= ':' NUMBER
   *
   * BONDS
   *
   *   bond_expression ::= bond_primitive
   *       | unary_operator bond_primitive
   *       | bond_expression bond_primitive
   *       | bond_expression binary_operator bond_primitive
   *
   *   bond_primitive ::= '-' | '=' | '#' | '$' | ':' | '/' | '\' | '~' | '@'
   *
   * CHAINS
   *
   *   ringbond ::= bond_expression? DIGIT | bond_expression? '%' DIGIT DIGIT
   *   
   *   branched_atom ::= atom ringbond* branch*
   *   
   *   branch ::= '(' chain ')' | '(' bond_expression chain ')' | '(' dot chain ')'
   *   
   *   chain ::= branched_atom | chain branched_atom | chain bond_expression branched_atom | chain dot branched_atom
   *   
   *   dot ::= '.'
   *
   * SMARTS STRINGS
   *
   *   smarts ::= chain terminator
   *   
   *   terminator ::= SPACE TAB | LINEFEED | CARRIAGE_RETURN | END_OF_STRING
   * @endcode
   *
   * @section license License
   * @code
   * Copyright (c) 2012, Tim Vandermeersch
   * All rights reserved.
   *
   * Redistribution and use in source and binary forms, with or without
   * modification, are permitted provided that the following conditions are met:
   *     * Redistributions of source code must retain the above copyright
   *       notice, this list of conditions and the following disclaimer.
   *     * Redistributions in binary form must reproduce the above copyright
   *       notice, this list of conditions and the following disclaimer in the
   *       documentation and/or other materials provided with the distribution.
   *     * Neither the name of the <organization> nor the
   *       names of its contributors may be used to endorse or promote products
   *       derived from this software without specific prior written permission.
   *
   * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
   * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
   * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
   * DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
   * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
   * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
   * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
   * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
   * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
   * @endcode
   */


  /**
   * The elements.
   */
  enum Elements {
    H = 1, He, Li, Be, B, C, N, O, F, Ne, Na, Mg, Al, Si, P, S, Cl, Ar, K, Ca, Sc, Ti, V,
    Cr, Mn, Fe, Co, Ni, Cu, Zn, Ga, Ge, As, Se, Br, Kr, Rb, Sr, Y, Zr, Nb, Mo, Tc, Ru,
    Rh, Pd, Ag, Cd, In, Sn, Sb, Te, I, Xe, Cs, Ba, Hf, Ta, W, Re, Os, Ir, Pt, Au, Hg,
    Tl, Pb, Bi, Po, At, Rn, Fr, Ra, Rf, Db, Sg, Bh, Hs, Mt, Ds, Rg, La, Ce, Pr, Nd, Pm,
    Sm, Eu, Gd, Tb, Dy, Ho, Er, Tm, Yb, Lu, Ac, Th, Pa, U, Np, Pu, Am, Cm, Bk, Cf, Es,
    Fm, Md, No, Lr
  };

  /**
   * Chirality classes.
   */
  enum Chirality {
    NotChiral = 0,
    AntiClockwise, Clockwise,
    TH1, TH2,
    AL1, AL2,
    SP1, SP2, SP3,
    TB1, TB2, TB3, TB4, TB5, TB6, TB7, TB8, TB9, TB10,
    TB11, TB12, TB13, TB14, TB15, TB16, TB17, TB18, TB19, TB20,
    OH1, OH2, OH3, OH4, OH5, OH6, OH7, OH8, OH9, OH10,
    OH11, OH12, OH13, OH14, OH15, OH16, OH17, OH18, OH19, OH20,
    OH21, OH22, OH23, OH24, OH25, OH26, OH27, OH28, OH29, OH30
  };

  /**
   * Error codes for Exception exceptions.
   */
  enum ErrorCode {
    ////////////////////////////////////////
    //
    // SyntaxError
    //
    ////////////////////////////////////////
    /**
     * Example: "[C"
     */
    NoClosingAtomBracket = 1,
    /**
     * Examples: "[]", "[13]", "[+]"
     */
    NoSymbolInBracketAtom = 2,
    /**
     * Examples: "[C@@T]", "[C@@A]", "[C@@TB]", "[C@@TH99]", "[C@@OH0]"
     */
    InvalidChirality = 3,
    /**
     * Example: "[C:]"
     */
    NoAtomClass = 4,
    /**
     * Example: "CC(C"
     */
    UnmatchedBranchOpening = 5,
    /**
     * Example: "CC)C"
     */
    UnmatchedBranchClosing = 6,
    /**
     * Examples: "Q", "!", "&", "Mm", "r"
     */
    InvalidAtomExpr = 7,
    /**
     * Example: "[13C$]", "[NP]"
     */
    TrailingCharInBracketAtom = 8,
    /**
     * Example: ".C"
     */
    LeadingDot = 9,
    /**
     * Example: "C."
     */
    TrailingDot = 10,
    /**
     * Example: "C%123CC%123", "C%CC%"
     */
    InvalidRingBondNumber = 11,
    ////////////////////////////////////////
    //
    // SemanticsError
    //
    ////////////////////////////////////////
    /**
     * Example: "[HH1]"
     */
    HydrogenHydrogenCount = 16,
    /**
     * Example: "C1CC"
     */
    UnmatchedRingBond = 32,
    /**
     * Example: "C-1CCCC=1"
     */
    ConflictingRingBonds = 64,
    /**
     * Example: "C12CCCCC12", "C11"
     */
    InvalidRingBond = 128,
    /**
     * Example: "C[C@H](F)(Cl)(Br)I", "O[C@]"
     */
    InvalidChiralValence = 256,
    /**
     * Example: "N[C@H2](F)I"
     */
    InvalidChiralHydrogenCount = 512
  };

  class Exception
  {
    public:
      /**
       * Exception type.
       */
      enum Type { SyntaxError, SemanticsError };

      /**
       * Constructor.
       *
       * @param type The exception type.
       * @param errorCode The numeric error code.
       * @param what Details of the exception.
       * @param pos The position in the SMILES/SMARTS string.
       * @param length The length of the error in the SMILES/SMARTS string.
       */
      Exception(Type type, ErrorCode errorCode, const std::string &what, std::size_t pos, std::size_t length = 1)
          : m_type(type), m_errorCode(errorCode), m_what(what), m_pos(pos), m_length(length)
      {
      }

      /**
       * The type of the exception.
       */
      Type type() const
      {
        return m_type;
      }

      /**
       * The ErrorCode for the exception.
       */
      ErrorCode errorCode() const
      {
        return m_errorCode;
      }

      /**
       * Details concerning the exception.
       */
      const std::string& what() const
      {
        return m_what;
      }

      /**
       * The position in the specified string where the error starts.
       */
      std::size_t pos() const
      {
        return m_pos;
      }

      /**
       * The length of the error.
       */
      std::size_t length() const
      {
        return m_length;
      }

    private:
      Type m_type;
      ErrorCode m_errorCode;
      std::string m_what;
      std::size_t m_pos;
      std::size_t m_length;
  };

  struct CallbackBase
  {
    void clear() {}
    void addAtom(int element, bool aromatic, int isotope, int hCount, int charge, int atomClass) {}
    void addBond(int source, int target, int order, bool isUp, bool isDown) {}
    void setChiral(int index, Chirality chirality, const std::vector<int> &chiralNbrs) {}
  };

  /**
   * Example Callback implementation to print the parsed results.
   */
  struct PrintCallback : public CallbackBase
  {
    /**
     * The clear() method is invoked when Parser::parse() is called and should
     * be used to initialize the callback function object to receive events for
     * a new SMILES/SMARTS.
     */
    void clear()
    {
    }

    /**
     * The addAtom() method is invoked when an atom is completly parsed.
     */
    void addAtom(int element, bool aromatic, int isotope, int hCount, int charge, int atomClass)
    {
      std::cout << "addAtom:" << std::endl
        << "    element: " << element << std::endl
        << "    aromatic: " << aromatic << std::endl
        << "    isotope: " << isotope << std::endl
        << "    hCount: " << hCount << std::endl
        << "    charge: " << charge << std::endl
        << "    atomClass: " << atomClass << std::endl;
    }

    /**
     * The addBond() method is invoked once both of the bond's atoms have been
     * added (by calling addAtom()). Therefore are the bond indices always
     * valid.
     */
    void addBond(int source, int target, int order, bool isUp, bool isDown)
    {
      std::cout << "addBond:" << std::endl
        << "    source: " << source << std::endl
        << "    target: " << target << std::endl
        << "    order: " << order << std::endl
        << "    siUp: " << isUp << std::endl
        << "    isDown: " << isDown << std::endl;
    }
    
    /**
     * The setChiral() method is invoked at the end of parsing for each atom
     * that has a chirality specified.
     */
    void setChiral(int index, Chirality chirality, const std::vector<int> &chiralNbrs)
    {
    }
  };

  /**
   * Value for implicit hydrogen in chiral nbrs.
   */
  inline int implicitHydrogen()
  {
    return std::numeric_limits<int>::max();
  }


  /**
   * @class Parser smiley.h <smiley.h>
   *
   * This is the main class for the Smiley SMILES/SMARTS parser.
   *
   * @code
   * #include "smiley.h"
   * 
   * using namespace Smiley;
   * 
   * int main(int argc, char **argv)
   * {
   *   PrintCallback callback;
   *   Parser<PrintCallback> parser(callback);
   * 
   *   try {
   *     parser.parse(argv[1]);
   *   } catch (Exception &e) {
   *     if (e.type() == Exception::SyntaxError)
   *       std::cerr << "Syntax";
   *     else
   *       std::cerr << "Semantics";
   *     std::cerr << "Error: " << e.what() << "." << std::endl;
   *     std::cerr << argv[1] << std::endl;
   *     for (std::size_t i = 0; i < e.pos(); ++i)
   *       std::cerr << " ";
   *     for (std::size_t i = 0; i < e.length(); ++i)
   *       std::cerr << "^";
   *     std::cerr << std::endl;
   *   }
   * }
   * @endcode
   */
  template<typename Callback>
  class Parser
  {
    private:
      /**
       * Internal structure to hold ring bond information.
       */
      struct RingBondInfo
      {
        RingBondInfo() : number(-1), order(-1), isUp(false), isDown(false),
            pos(std::string::npos)
        {
        }

        RingBondInfo(int number_, int order_, bool isUp_, bool isDown_,
            bool isExplicit_, std::size_t pos_) : number(number_),
            order(order_), isUp(isUp_), isDown(isDown_),
            isExplicit(isExplicit_), pos(pos_)
        {
        }

        std::size_t pos;
        int number;
        int order;
        bool isUp;
        bool isDown;
        bool isExplicit;
      };

      /**
       * Internal structure to hold branch information.
       */
      struct BranchInfo
      {
        BranchInfo() : pos(std::string::npos), index(-1)
        {
        }

        BranchInfo(int index_, std::size_t pos_) : pos(pos_), index(index_)
        {
        }

        std::size_t pos;
        int index;
      };

      /**
       * Internal structure to hold chiral information
       */
      struct ChiralInfo
      {
        ChiralInfo() : pos(std::string::npos), chiral(NotChiral)
        {
        }

        ChiralInfo(Chirality chirality, std::size_t pos_) : pos(pos_), chiral(chirality)
        {
        }

        std::size_t pos;
        std::vector<int> nbrs;
        Chirality chiral;
      };

      /**
       * find a matching bracket with the open bracket at position @p pos.
       *
       * @param open Open bracket symbol (e.g. "[" or "(").
       * @param close Close bracket symbol (e.g. "]" or ")").
       * @param pos Position of the opening bracket.
       * @return The position of the matching closing bracket. Throws SyntaxError
       *         exception if there is no matching bracket.
       *
       * Examples:
       * @verbatim
       * findMatchingBracket("[", "]", 2) -> 6, 15
       *
       * str: "CC[13C]CC"     str: "CC[N$(CC[NH3+])]CC"
       *         ^   ^                ^            ^
       *         2   6                2            15
       * @endverbatim
       */
      std::size_t findMatchingBracket(const std::string &open, const std::string &close, std::size_t pos = 0)
      {
        std::size_t indent = 1;
        while (indent) {
          std::size_t open_pos = m_str.find(open, pos + 1);
          std::size_t close_pos = m_str.find(close, pos + 1);
          // syntax error when there is no closing bracket
          if (close_pos == std::string::npos)
            throw Exception(Exception::SyntaxError, NoClosingAtomBracket, 
                "Could not find matching bracket", pos, m_str.size() - pos);
          if (open_pos > close_pos) {
            --indent;
            pos = close_pos;
          } else {
            ++indent;
            pos = open_pos;
          }
        }
        return pos;
      }

      /**
       * Check the character at position @p m_pos. If the next character matches
       * @p chr, position @p m_pos is incremeneted.
       *
       * @param chr The character to check for.
       * @return True if the character matches.
       */
      bool checkNextChar(char chr)
      {
        if (m_pos + 1 >= m_str.size())
          return false;
        bool match = m_str[m_pos + 1] == chr;
        if (match)
          ++m_pos;
        return match;
      }

      /**
       * Add a bond by calling the callback's addBond() member function.
       * 
       * @param source Bond source atom index (0...N)
       * @param target Bond target atom index (0...N)
       * @param order The bond order.
       */
      void addBond(int source, int target, int order, bool isUp, bool isDown)
      {
        for (std::size_t i = 0; i < m_chiralInfo[source].nbrs.size(); i++) {
          int nbr = m_chiralInfo[source].nbrs[i];
          if (nbr == target) {
            if (m_exceptions & InvalidRingBond)
              throw Exception(Exception::SemanticsError, InvalidRingBond,
                  "Parallel ring bond", 0, 0);
            else
              return; // ignore parallel bond
          }
        }
        if (source == target) {
          if (m_exceptions & InvalidRingBond)
            throw Exception(Exception::SemanticsError, InvalidRingBond,
                "Self-loop ring bond", 0, 0);
          else
            return; // ignore self-loop
        }

        m_callback.addBond(source, target, order, isUp, isDown);

        m_chiralInfo[source].nbrs.push_back(target);
        if (m_chiralInfo[target].nbrs.size() && m_chiralInfo[target].nbrs.front() == implicitHydrogen())
          m_chiralInfo[target].nbrs.insert(m_chiralInfo[target].nbrs.begin(), source);
        else
          m_chiralInfo[target].nbrs.push_back(source);
      }

      /**
       * Add an atom by calling the callback's addAtom() member function. If
       * there is a previous atom (i.e. @p m_prev != -1) then the bond between
       * the new atom and the previous atom will by added by calling addBond()
       * after adding the new atom.
       */
      void addAtom(int element, bool aromatic = false, int isotope = -1, int hCount = -1, int charge = 0, int atomClass = 0)
      {
        // check for [HH1]
        if (element == H && hCount != 0)
          throw Exception(Exception::SemanticsError, HydrogenHydrogenCount,
              "Hydrogen atoms can not have a hydrogen count", 0, 0);

        m_callback.addAtom(element, aromatic, isotope, hCount, charge, atomClass);

        if (m_prev != -1)
          addBond(m_prev, m_index, m_bondOrder, m_isUp, m_isDown);

        m_prev = m_index;
        ++m_index;
        m_chiralInfo.push_back(ChiralInfo());
      }

      /**
       * Reset the bond information.
       */
      void resetBondInfo()
      {
        // reset bond info
        m_bondOrder = 1;
        m_isUp = false;
        m_isDown = false;
        m_explicitBond = false;
      }

      /**
       * @code
       * isotope ::= NUMBER
       * @endcode
       */
      void parseIsotope()
      {
        m_isotope = 0;
        bool found_isotope = false;
        while (std::isdigit(m_str[m_pos])) {
          m_isotope *= 10;
          m_isotope += m_str[m_pos] - '0';
          ++m_pos;
          found_isotope = true;
        }
        if (!found_isotope)
          m_isotope = -1;
      }

      /**
       * @code
       * symbol ::= element_symbols | aromatic_symbols | '*'
       * element_symbols ::= 'H' | 'He' | 'Li' | ... | 'No' | 'Lr'
       * aromatic_symbols ::= 'c' | 'n' | 'o' | 'p' | 's' | 'se' | 'as'
       * @endcode
       */
      void parseSymbol()
      {
        switch (m_str[m_pos]) {
          case 'H':
            if (checkNextChar('e'))
              m_element = He;
            else if (checkNextChar('f'))
              m_element = Hf;
            else if (checkNextChar('g'))
              m_element = Hg;
            else if (checkNextChar('s'))
              m_element = Hs;
            else if (checkNextChar('o'))
              m_element = Ho;
            else
              m_element = H;
            break;
          case 'L':
            if (checkNextChar('i'))
              m_element = Li;
            else if (checkNextChar('a'))
              m_element = La;
            else if (checkNextChar('u'))
              m_element = Lu;
            else if (checkNextChar('r'))
              m_element = Lr;
            break;
          case 'B':
            if (checkNextChar('e'))
              m_element = Be;
            else if (checkNextChar('r'))
              m_element = Br;
            else if (checkNextChar('a'))
              m_element = Ba;
            else if (checkNextChar('i'))
              m_element = Bi;
            else if (checkNextChar('h'))
              m_element = Bh;
            else if (checkNextChar('k'))
              m_element = Bk;
            else
              m_element = B;
            break;
          case 'C':
            if (checkNextChar('l'))
              m_element = Cl;
            else if (checkNextChar('a'))
              m_element = Ca;
            else if (checkNextChar('r'))
              m_element = Cr;
            else if (checkNextChar('o'))
              m_element = Co;
            else if (checkNextChar('u'))
              m_element = Cu;
            else if (checkNextChar('d'))
              m_element = Cd;
            else if (checkNextChar('s'))
              m_element = Cs;
            else if (checkNextChar('e'))
              m_element = Ce;
            else if (checkNextChar('m'))
              m_element = Cm;
            else if (checkNextChar('f'))
              m_element = Cf;
            else
              m_element = C;
            break;
          case 'N':
            if (checkNextChar('e'))
              m_element = Ne;
            else if (checkNextChar('a'))
              m_element = Na;
            else if (checkNextChar('i'))
              m_element = Ni;
            else if (checkNextChar('b'))
              m_element = Nb;
            else if (checkNextChar('d'))
              m_element = Nd;
            else if (checkNextChar('p'))
              m_element = Np;
            else if (checkNextChar('o'))
              m_element = No;
            else
              m_element = N;
            break;
          case 'O':
            if (checkNextChar('s'))
              m_element = Os;
            else
              m_element = O;
            break;
          case 'F':
            if (checkNextChar('e'))
              m_element = Fe;
            else if (checkNextChar('r'))
              m_element = Fr;
            else if (checkNextChar('m'))
              m_element = Fm;
            else
              m_element = F;
            break;
          case 'M':
            if (checkNextChar('n'))
              m_element = Mn;
            else if (checkNextChar('o'))
              m_element = Mo;
            else if (checkNextChar('t'))
              m_element = Mt;
            else if (checkNextChar('d'))
              m_element = Md;
            else if (checkNextChar('g'))
              m_element = Mg;
            break;
          case 'A':
            if (checkNextChar('l'))
              m_element = Al;
            else if (checkNextChar('r'))
              m_element = Ar;
            else if (checkNextChar('s'))
              m_element = As;
            else if (checkNextChar('g'))
              m_element = Ag;
            else if (checkNextChar('u'))
              m_element = Au;
            else if (checkNextChar('t'))
              m_element = At;
            else if (checkNextChar('c'))
              m_element = Ac;
            else if (checkNextChar('m'))
              m_element = Am;
            break;
          case 'S':
            if (checkNextChar('i'))
              m_element = Si;
            else if (checkNextChar('c'))
              m_element = Sc;
            else if (checkNextChar('e'))
              m_element = Se;
            else if (checkNextChar('r'))
              m_element = Sr;
            else if (checkNextChar('n'))
              m_element = Sn;
            else if (checkNextChar('b'))
              m_element = Sb;
            else if (checkNextChar('g'))
              m_element = Sg;
            else if (checkNextChar('m'))
              m_element = Sm;
            else
              m_element = S;
            break;
          case 'P':
            if (checkNextChar('d'))
              m_element = Pd;
            else if (checkNextChar('t'))
              m_element = Pt;
            else if (checkNextChar('b'))
              m_element = Pb;
            else if (checkNextChar('o'))
              m_element = Po;
            else if (checkNextChar('r'))
              m_element = Pr;
            else if (checkNextChar('m'))
              m_element = Pm;
            else if (checkNextChar('a'))
              m_element = Pa;
            else if (checkNextChar('u'))
              m_element = Pu;
            else
              m_element = P;
            break;
          case 'K':
            if (checkNextChar('r'))
              m_element = Kr;
            else
              m_element = K;
            break;
          case 'T':
            if (checkNextChar('i'))
              m_element = Ti;
            else if (checkNextChar('c'))
              m_element = Tc;
            else if (checkNextChar('e'))
              m_element = Te;
            else if (checkNextChar('a'))
              m_element = Ta;
            else if (checkNextChar('l'))
              m_element = Tl;
            else if (checkNextChar('b'))
              m_element = Tb;
            else if (checkNextChar('m'))
              m_element = Tm;
            else if (checkNextChar('h'))
              m_element = Th;
            break;
          case 'V':
            m_element = V;
            break;
          case 'Z':
            if (checkNextChar('n'))
              m_element = Zn;
            else if (checkNextChar('r'))
              m_element = Zr;
            break;
          case 'G':
            if (checkNextChar('a'))
              m_element = Ga;
            else if (checkNextChar('e'))
              m_element = Ge;
            else if (checkNextChar('d'))
              m_element = Gd;
            break;
          case 'R':
            if (checkNextChar('b'))
              m_element = Rb;
            else if (checkNextChar('u'))
              m_element = Ru;
            else if (checkNextChar('h'))
              m_element = Rh;
            else if (checkNextChar('e'))
              m_element = Re;
            else if (checkNextChar('n'))
              m_element = Rn;
            else if (checkNextChar('a'))
              m_element = Ra;
            else if (checkNextChar('f'))
              m_element = Rf;
            else if (checkNextChar('g'))
              m_element = Rg;
            break;
          case 'Y':
            if (checkNextChar('b'))
              m_element = Yb;
            else
              m_element = Y;
            break;
          case 'I':
            if (checkNextChar('n'))
              m_element = In;
            else if (checkNextChar('r'))
              m_element = Ir;
            else
              m_element = I;
            break;
          case 'X':
            if (checkNextChar('e'))
              m_element = Xe;
            break;
          case 'W':
            m_element = W;
            break;
          case 'D':
            if (checkNextChar('b'))
              m_element = Db;
            else if (checkNextChar('s'))
              m_element = Ds;
            else if (checkNextChar('y'))
              m_element = Dy;
            break;
          case 'E':
            if (checkNextChar('u'))
              m_element = Eu;
            else if (checkNextChar('r'))
              m_element = Er;
            else if (checkNextChar('s'))
              m_element = Es;
            break;
          case 'U':
            m_element = U;
            break;
          case 'c':
            m_element = C;
            m_aromatic = true;
            break;
          case 'n':
            m_element = N;
            m_aromatic = true;
            break;
          case 'o':
            m_element = O;
            m_aromatic = true;
            break;
          case 'p':
            m_element = P;
            m_aromatic = true;
            break;
          case 's':
            if (checkNextChar('e'))
              m_element = Se;
            else
              m_element = S;
            m_aromatic = true;
            break;
          case 'a':
            if (checkNextChar('s')) {
              m_element = As;
              m_aromatic = true;          
            }
            break;
        }

        if (!m_element && m_str[m_pos] != '*')
          throw Exception(Exception::SyntaxError, NoSymbolInBracketAtom,
              "Bracket atom expression does not contain element symbol", m_pos, 1);
        else
          ++m_pos;
      }

      /**
       * @code
       * chiral ::= '@@' | '@@@@' | '@@TH1' | '@@TH2' | '@@AL1' | '@@AL2'
       *     | '@@SP1' | '@@SP2' | '@@SP3'
       *     | '@@TB1' | '@@TB2' | ... | '@@TB20'
       *     | '@@OH1' | '@@OH2' | ... | '@@OH30'
       * @endcode
       */
      void parseChiral()
      {
        m_chiralInfo.back().pos = m_pos;

        if (m_str[m_pos] != '@')
          return;

        // @@
        if (checkNextChar('@')) {
          m_chiral = Clockwise;
          ++m_pos;
          return;
        }
        if (checkNextChar('T')) {
          // @TH1 & @TH2
          if (checkNextChar('H')) {
            if (checkNextChar('1'))
              m_chiral = TH1;
            else if (checkNextChar('2'))
              m_chiral = TH2;
            else
              throw Exception(Exception::SyntaxError, InvalidChirality,
                  "Invalid chiral specified, expected 1 or 2", m_pos + 1, 1);
            ++m_pos;
            return;
          }
          // @TB1 ... @TB20
          if (checkNextChar('B')) {
            std::size_t pos = m_pos;
            int tb = 0;
            if (std::isdigit(m_str[m_pos + 1])) {
              tb = m_str[m_pos + 1] - '0';
              ++m_pos;
            }
            if (std::isdigit(m_str[m_pos + 1])) {
              tb *= 10;
              tb += m_str[m_pos + 1] - '0';
              ++m_pos;
            }
            if (tb < 1 || tb > 20)
              throw Exception(Exception::SyntaxError, InvalidChirality,
                  "Invalid chiral class specified, expected 1-20", pos + 1, m_pos == pos ? 1 : m_pos - pos);

            m_chiral = static_cast<Chirality>(tb + TB1 - 1);
            ++m_pos;
            return;
          }

          throw Exception(Exception::SyntaxError, InvalidChirality,
              "Invalid chiral specifier, expected H or B", m_pos + 1, 1);
        }
        // @AL1 & @AL2
        if (checkNextChar('A')) {
          if (!checkNextChar('L'))
            throw Exception(Exception::SyntaxError, InvalidChirality,
                "Invalid chiral specifier, expected L", m_pos + 1, 1);
          if (checkNextChar('1'))
            m_chiral = AL1;
          else if (checkNextChar('2'))
            m_chiral = AL2;
          else
            throw Exception(Exception::SyntaxError,  InvalidChirality,
                "Invalid chiral specified, expected 1 or 2", m_pos + 1, 1);
          ++m_pos;
          return;
        }
        // @SP1, @SP2 & @SP3
        if (checkNextChar('S')) {
          if (!checkNextChar('P'))
            throw Exception(Exception::SyntaxError, InvalidChirality,
                "Invalid chiral specifier, expected P", m_pos + 1, 1);
          if (checkNextChar('1'))
            m_chiral = SP1;
          else if (checkNextChar('2'))
            m_chiral = SP2;
          else if (checkNextChar('3'))
            m_chiral = SP3;
          else
            throw Exception(Exception::SyntaxError, InvalidChirality,
                "Invalid chiral specified, expected 1, 2 or 3", m_pos + 1, 1);
          ++m_pos;
          return;
        }
        // @OH1 ... @OH30
        if (checkNextChar('O')) {
          if (!checkNextChar('H'))
            throw Exception(Exception::SyntaxError, InvalidChirality,
                "Invalid chiral specifier, expected H", m_pos + 1, 1);

          std::size_t pos = m_pos;
          int oh = 0;
          if (std::isdigit(m_str[m_pos + 1])) {
            oh = m_str[m_pos + 1] - '0';
            ++m_pos;
          }
          if (std::isdigit(m_str[m_pos + 1])) {
            oh *= 10;
            oh += m_str[m_pos + 1] - '0';
            ++m_pos;
          }
          if (oh < 1 || oh > 30)
            throw Exception(Exception::SyntaxError, InvalidChirality,
                "Invalid chiral class specified, expected 1-30", pos + 1, m_pos == pos ? 1 : m_pos - pos);

          m_chiral = static_cast<Chirality>(oh + OH1 - 1);
          ++m_pos;
          return;
        }

        m_chiral = AntiClockwise;
        ++m_pos;
      }

      /**
       * @code
       * hcount ::= 'H' | 'H' DIGIT
       * @endcode
       */
      void parseHydrogenCount()
      {
        if (DEBUG)
          std::cout << "parseHydrogenCount(" << m_str.substr(m_pos) << ")" << std::endl;

        // [C] = [CH0]
        m_hCount = 0;
        if (m_str[m_pos] != 'H')
          return;
        ++m_pos;
        if (std::isdigit(m_str[m_pos])) {
          m_hCount = m_str[m_pos] - '0';
          ++m_pos;
          return;
        }
        m_hCount = 1;
      }

      /**
       * @code
       * charge ::= '-' | '+' | '-' DIGIT | '+' DIGIT | '--' | '++'
       * @endcode
       *
       * Deprecated: '--' and '++'
       */
      void parseCharge()
      {
        if (DEBUG)
          std::cout << "parseCharge(" << m_str.substr(m_pos) << ")" << std::endl;

        if (m_str[m_pos] == '-') {
          if (checkNextChar('-')) {
            m_charge = -2;
            ++m_pos;
            return;
          }
          if (std::isdigit(m_str[m_pos + 1])) {
            m_charge = - (m_str[m_pos + 1] - '0');
            m_pos += 2;
            return;
          }
          ++m_pos;
          m_charge = -1;
        } else if (m_str[m_pos] == '+') {
          if (checkNextChar('+')) {
            m_charge = 2;
            ++m_pos;
            return;        
          }
          if (std::isdigit(m_str[m_pos + 1])) {
            m_charge = m_str[m_pos + 1] - '0';
            m_pos += 2;
            return;
          }
          ++m_pos;
          m_charge = 1;
        }
      }

      /**
       * @code
       * class ::= ':' NUMBER
       * @endcode
       */
      void parseClass()
      {
        if (DEBUG)
          std::cout << "parseClass(" << m_str.substr(m_pos) << ")" << std::endl;

        if (m_str[m_pos] != ':')
          return;
        bool found_number = false;
        while (std::isdigit(m_str[m_pos + 1])) {
          m_class *= 10;
          m_class += m_str[m_pos + 1] - '0';
          ++m_pos;
          found_number = true;
        }
        ++m_pos;
        if (!found_number)
          throw Exception(Exception::SyntaxError, NoAtomClass,
              "No atom class, expected number", m_pos + 1, 1);
      }

      /**
       * @code
       * bracket_atom ::= '[' isotope? symbol chiral? hcount? charge? class? ']'
       * @endcode
       */
      void parseBracketAtom()
      {
        if (DEBUG)
          std::cout << "parseBracketAtom(" << m_str.substr(m_pos) << ")" << std::endl;

        std::size_t close = findMatchingBracket("[", "]", m_pos);
        ++m_pos;

        parseIsotope();
        //std::cout << "isotope: " << m_isotope << std::endl;
        parseSymbol();
        //std::cout << "element: " << m_element << std::endl;
        parseChiral();
        //std::cout << "chrial: " << m_chiral << std::endl;
        parseHydrogenCount();
        //std::cout << "hCount: " << m_hCount << std::endl;
        parseCharge();
        //std::cout << "charge: " << m_charge << std::endl;
        parseClass();
        //std::cout << "class: " << m_class << std::endl;

        m_chiralInfo.back().chiral = m_chiral;
        if (m_hCount >= 1)
          m_chiralInfo.back().nbrs.push_back(implicitHydrogen());
        if (m_hCount > 1 && m_chiral && m_exceptions & InvalidChiralHydrogenCount) {
          throw Exception(Exception::SemanticsError, InvalidChiralHydrogenCount,
              "Chiral atoms can only have one hydrogen", m_chiralInfo.back().pos, 1);
        } 

        //std::cout << m_str[m_pos] << std::endl;
        if (m_str[m_pos] != ']')
          throw Exception(Exception::SyntaxError, TrailingCharInBracketAtom,
              "Bracket atom expression contains invalid trailing characters", m_pos, close - m_pos);

        addAtom(m_element, m_aromatic, m_isotope, m_hCount, m_charge, m_class);
      }

      /**
       * @code
       * aliphatic_organic ::= 'B' | 'C' | 'N' | 'O' | 'P' | 'S' | 'F' | 'Cl' | 'Br' | 'I'
       * aromatic_organic ::= 'b' | 'c' | 'n' | 'o' | 'p' | 's'
       * @endcode
       */
      bool parseOrganicSubsetAtom()
      {
        if (DEBUG)
          std::cout << "parseOrganicSubsetAtom(" << m_str.substr(m_pos) << ")" << std::endl;

        switch (m_str[m_pos]) {
          case 'B':
            if (checkNextChar('r'))
              addAtom(Br);
            else
              addAtom(B);
            ++m_pos;
            return true;
          case 'C':
            if (checkNextChar('l'))
              addAtom(Cl);
            else
              addAtom(C);
            ++m_pos;
            return true;
          case 'N':
            addAtom(N);
            ++m_pos;
            return true;
          case 'O':
            addAtom(O);
            ++m_pos;
            return true;
          case 'S':
            addAtom(S);
            ++m_pos;
            return true;
          case 'P':
            addAtom(P);
            ++m_pos;
            return true;
          case 'F':
            addAtom(F);
            ++m_pos;
            return true;
          case 'I':
            addAtom(I);
            ++m_pos;
            return true;
          case 'b':
            addAtom(B, true);
            ++m_pos;
            return true;
          case 'c':
            addAtom(C, true);
            ++m_pos;
            return true;
          case 'n':
            addAtom(N, true);
            ++m_pos;
            return true;
          case 'o':
            addAtom(O, true);
            ++m_pos;
            return true;
          case 's':
            addAtom(S, true);
            ++m_pos;
            return true;
          case 'p':
            addAtom(P, true);
            ++m_pos;
            return true;
        }

        return false;
      }

      /**
       * @code
       * atom ::= bracket_atom | aliphatic_organic | aromatic_organic | '*'
       * @endcode
       */
      bool parseAtom()
      {
        if (DEBUG)
          std::cout << "parseAtom(" << m_str.substr(m_pos) << ")" << std::endl;

        m_element = 0;
        m_isotope = -1;
        m_charge = 0;
        m_chiral = NotChiral;
        m_hCount = -1;
        m_class = 0;
        m_aromatic = false;

        switch (m_str[m_pos]) {
          case '[':
            parseBracketAtom();
            ++m_pos;
            return true;
          case '*':
            addAtom(0);
            ++m_pos;
            return true;
        }

        return parseOrganicSubsetAtom();
      }

      /**
       * @code
       * bond ::= '-' | '=' | '#' | '$' | ':' | '/' | '\'
       * @endcode
       */
      void parseBond()
      {
        if (m_pos >= m_str.size())
          return;
        if (DEBUG)
          std::cout << "parseBond(" << m_str.substr(m_pos) << ")" << std::endl;

        if (m_str[m_pos] == '-') {
          m_bondOrder = 1;
          m_explicitBond = true;
          ++m_pos;
        } else if (m_str[m_pos] == '=') {
          m_bondOrder = 2;
          m_explicitBond = true;
          ++m_pos;
        } else if (m_str[m_pos] == '#') {
          m_bondOrder = 3;
          m_explicitBond = true;
          ++m_pos;
        } else if (m_str[m_pos] == '$') {
          m_bondOrder = 4;
          m_explicitBond = true;
          ++m_pos;
        } else if (m_str[m_pos] == ':') {
          m_bondOrder = 5;
          m_explicitBond = true;
          ++m_pos;
        } else if (m_str[m_pos] == '/') {
          m_isUp = true;
          m_explicitBond = true;
          ++m_pos;
        } else if (m_str[m_pos] == '\\') {
          m_isDown = true;
          m_explicitBond = true;
          ++m_pos;
        }
      }

      /**
       * Used for debugging.
       */
      void printRingBonds()
      {
        for (typename std::map<int, std::vector<RingBondInfo> >::iterator i = m_ringBonds.begin(); i != m_ringBonds.end(); ++i) {
          std::cout << "    RingBond: index = " << i->first << std::endl;
          for (std::size_t j = 0; j < i->second.size(); ++j)
            std::cout << "        " << i->second[j].number << std::endl; 
        }
      }

      /**
       * Helper function for parseRingBond().
       */
      void processRingBond(int rnum, std::size_t pos)
      {
        //std::cout << "BEFORE processing " << rnum << std::endl; printRingBonds();
        typename std::map<int, std::vector<RingBondInfo> >::iterator ringBond = m_ringBonds.begin();
        std::size_t j;
        bool found = false;
        for (; ringBond != m_ringBonds.end(); ++ringBond) {
          for (j = 0; j < ringBond->second.size(); ++j)
            if (ringBond->second[j].number == rnum) {
              found = true;
              break;
            }
          if (found)
            break;
        }
        if (ringBond != m_ringBonds.end()) {
          if (ringBond->second[j].isExplicit) {
            if (m_explicitBond && m_exceptions & ConflictingRingBonds) {
              if (ringBond->second[j].order != m_bondOrder || 
                  ringBond->second[j].isUp != m_isUp ||
                  ringBond->second[j].isDown != m_isDown)
                throw Exception(Exception::SemanticsError, ConflictingRingBonds,
                    "Conflicing ring bonds", pos, 1);
            }
            addBond(ringBond->first, m_prev, ringBond->second[j].order, ringBond->second[j].isUp, ringBond->second[j].isDown);
          } else
            addBond(ringBond->first, m_prev, m_bondOrder, m_isUp, m_isDown);
          ringBond->second.erase(ringBond->second.begin() + j);
          if (ringBond->second.empty())
            m_ringBonds.erase(ringBond);
        } else
          m_ringBonds[m_prev].push_back(RingBondInfo(rnum, m_bondOrder, m_isUp, m_isDown, m_explicitBond, pos));
        //std::cout << "AFTER processing " << rnum << std::endl; printRingBonds();
        resetBondInfo();
      }

      /**
       * @code
       * ringbond ::= bond? DIGIT | bond? '%' DIGIT DIGIT
       * @endcode
       */
      void parseRingBond()
      {
        if (DEBUG)
          std::cout << "parseRingBond(" << m_str.substr(m_pos) << ")" << std::endl;

        parseBond();
        if (std::isdigit(m_str[m_pos])) {
          processRingBond(m_str[m_pos] - '0', m_pos);
          ++m_pos;
        } else if (m_str[m_pos] == '%') {
          std::size_t bond_pos = m_pos - 1;
          if (m_pos + 2 >= m_str.size())
            throw Exception(Exception::SyntaxError, InvalidRingBondNumber,
                "Invalid ring bond, expected number", m_pos + 1, 2);
          if (!std::isdigit(m_str[m_pos + 1]) || !std::isdigit(m_str[m_pos + 2]))
            throw Exception(Exception::SyntaxError, InvalidRingBondNumber,
                "Expected ring bond number", m_pos + 1, 2);
          int r = 10 * (m_str[m_pos + 1] - '0') + m_str[m_pos + 2] - '0';
          processRingBond(r, bond_pos);
          m_pos += 3;
        }
      }

      /**
       * @code
       * dot ::= '.'
       * terminator ::= SPACE TAB | LINEFEED | CARRIAGE_RETURN | END_OF_STRING
       * branched_atom ::= atom ringbond* branch*
       * branch ::= '(' chain ')' | '(' bond chain ')' | '(' dot chain ')'
       * chain ::= branched_atom | chain branched_atom | chain bond branched_atom | chain dot branched_atom
       * smiles ::= chain terminator
       * @endcode
       */
      void parseChain()
      {
        if (DEBUG)
          std::cout << "parseChain(" << m_str.substr(m_pos) << ")" << std::endl;

        while (true) {
          // check for dot ::= '.'?
          if (m_str[m_pos] == '.') {
            if (m_index == 0)
              throw Exception(Exception::SyntaxError, LeadingDot,
                  "Found dot '.' at begining of pattern", 0, 1);
            if (m_pos + 1 >= m_str.size())
              throw Exception(Exception::SyntaxError, TrailingDot,
                  "Found dor '.' at ending of pattern", m_pos - 1, 1);
            ++m_pos;
            m_prev = -1;
          }

          // check for closing branch ::= ')'*
          while (m_str[m_pos] == ')') {
            if (m_branches.size()) {
              if (DEBUG)
                std::cout << "    close branch: " << m_branches.back().index << " @pos " << m_pos << std::endl;
              m_prev = m_branches.back().index;
              m_branches.pop_back();
              ++m_pos;
            } else
              throw Exception(Exception::SyntaxError, UnmatchedBranchClosing,
                  "Unmatched branch closing", 0, m_pos + 1);
            if (m_pos >= m_str.size())
              break;
          }

          // if there is a previous atom, check for bond ::= bond?
          if (m_prev != -1)
            parseBond();
          if (m_pos >= m_str.size())
            break;

          // check atom ::= atom
          if (!parseAtom() && m_str[m_pos] != '(')
            throw Exception(Exception::SyntaxError, InvalidAtomExpr,
                "Could not parse atom expression", m_pos);
          else
            resetBondInfo();
          if (m_pos >= m_str.size())
            break;

          // check for ring bond ::= ring_bond*
          std::size_t pos = std::string::npos;
          while (pos != m_pos) {
            pos = m_pos;
            parseRingBond();
          }
          if (m_pos >= m_str.size())
            break;

          // check for branch opening ::= '('?
          pos = std::string::npos;
          while (pos != m_pos) {
            pos = m_pos;
            if (m_str[m_pos] == '(') {
              if (DEBUG)
                std::cout << "    open branch: " << m_prev << " @pos " << m_pos << std::endl;
              m_branches.push_back(BranchInfo(m_prev, m_pos));
              ++m_pos;
              parseChain();
            }
          }

          if (m_pos >= m_str.size())
            break;

          // check for termination
          bool terminate = false;
          switch (m_str[m_pos]) {
            case ' ':
            case '\t':
            case '\n':
            case '\r':
              terminate = true;
              break;
            default:
              break;
          }
          if (terminate)
            break;
        }
      }

      /**
       * @return Expected valence.
       */
      int processAlleneStereochemistry(ChiralInfo &chiralInfo)
      {
        if (chiralInfo.nbrs.size() != 2)
          return -1;
        std::vector<int> &nbrs = chiralInfo.nbrs;
        int lft = nbrs[0];
        int rgt = nbrs[1];
        if (m_chiralInfo[lft].nbrs.size() != 3 || m_chiralInfo[rgt].nbrs.size() != 3)
          return -1;
        // NC(F)=[C@AL1]=C(Cl)I     nbrs[1] = [ 0 2 3 ]  <-----------+
        // ^  ^   ^        ^  ^     nbrs[3] = [ 1 4 ]    <- chiral   +-- neigbors
        // 0  2   3        5  6     nbrs[4] = [ 3 5 6 ]  <-----------+
        nbrs.clear();
        // copy [ 0 2 ] to nbrs[3]
        nbrs.insert(nbrs.end(), m_chiralInfo[lft].nbrs.begin(), m_chiralInfo[lft].nbrs.end() - 1);
        // copy [ 5 6 ] to nbrs[3]
        nbrs.insert(nbrs.end(), m_chiralInfo[rgt].nbrs.begin() + 1, m_chiralInfo[rgt].nbrs.end());
        if (chiralInfo.chiral == AntiClockwise)
          chiralInfo.chiral = AL1;
        else if (chiralInfo.chiral == Clockwise)
          chiralInfo.chiral = AL2;
        return 4;
      }

      /**
       * Handle stereochemistry by calling Callback::setChiral() for each chiral atom.
       */
      void processStereochemistry()
      {
        for (std::size_t i = 0; i < m_chiralInfo.size(); ++i) {
          if (!m_chiralInfo[i].chiral)
            continue;
          int expectedValence = -1;
          switch (m_chiralInfo[i].chiral) {
            case Clockwise:
              switch (m_chiralInfo[i].nbrs.size()) {
                case 2:
                  // Allene
                  expectedValence = processAlleneStereochemistry(m_chiralInfo[i]);
                  break;
                case 4:
                  // Tetrahedral
                  // Square Planar
                  expectedValence = 4;
                  // determining correct class requires knowledge of atom (TH vs. SP)
                  break;
                case 5:
                  m_chiralInfo[i].chiral = TB1;
                  expectedValence = 5;
                  break;
                case 6:
                  m_chiralInfo[i].chiral = OH1;
                  expectedValence = 6;
                  break;
              }
              break;
            case AntiClockwise:
              switch (m_chiralInfo[i].nbrs.size()) {
                case 2:
                  // Allene
                  expectedValence = processAlleneStereochemistry(m_chiralInfo[i]);
                  break;
                case 4:
                  expectedValence = 4;
                  break;
                case 5:
                  m_chiralInfo[i].chiral = TB2;
                  expectedValence = 5;
                  break;
                case 6:
                  m_chiralInfo[i].chiral = OH2;
                  expectedValence = 6;
                  break;
              }
              break;
            case TH1:
            case TH2:
            case SP1:
            case SP2:
            case SP3:
              expectedValence = 4;
              break;
            case AL1:
            case AL2:
              expectedValence = processAlleneStereochemistry(m_chiralInfo[i]);
              break;
            default:
              if (m_chiralInfo[i].chiral >= TB1 && m_chiralInfo[i].chiral <= TB20)
                expectedValence = 5;
              else if (m_chiralInfo[i].chiral >= OH1 && m_chiralInfo[i].chiral <= OH30)
                expectedValence = 6;
              break;
          }

          if (m_chiralInfo[i].nbrs.size() != expectedValence && m_exceptions & InvalidChiralValence)
            throw Exception(Exception::SemanticsError, InvalidChiralValence,
                "Invalid chiral valence", m_chiralInfo[i].pos, 1);

          m_callback.setChiral(i, m_chiralInfo[i].chiral, m_chiralInfo[i].nbrs);
        }
      }

    public:
      /**
       * Constructor.
       *
       * @param callback The callback function object to handle the information
       * resulting from parsing.
       */
      Parser(Callback &callback) : m_callback(callback), m_exceptions(~0)
      {
      }

      /**
       * Parse a SMILES/SMARTS string.
       *
       * @param str The SMILES/SMARTS string.
       */
      void parse(const std::string &str)
      {
        m_callback.clear();

        if (str.empty())
          return;

        m_str = str;
        m_pos = 0;

        m_index = 0;
        m_prev = -1;
        m_branches.clear();
        m_ringBonds.clear();
        m_chiralInfo.clear();
        m_chiralInfo.push_back(ChiralInfo());

        parseChain();

        if (m_branches.size())
          throw Exception(Exception::SyntaxError, UnmatchedBranchOpening,
              "Unmatched branch opening", m_branches.back().pos, m_str.size() - m_branches.back().pos);
        if (m_ringBonds.size() && m_exceptions & UnmatchedRingBond)
          throw Exception(Exception::SemanticsError, UnmatchedRingBond,
              "Unmatched ring bond", m_ringBonds.begin()->second[0].pos, 1);

        processStereochemistry();
      }

      /**
       * Disable the specified @p exceptions.
       *
       * @param exceptions Exception ErrorCode (may be OR-ed).
       */
      void disableExceptions(ErrorCode exceptions)
      {
        m_exceptions &= ~exceptions;
      }

      /**
       * Enable specified exceptions.
       *
       * @param exceptions Exceptions to enable (may be OR-ed).
       */
      void enableExceptions(ErrorCode exceptions)
      {
        m_exceptions |= exceptions;
      }

      /**
       * Disable all disableable exceptions.
       */
      void disableAllExceptions()
      {
        m_exceptions = 0;
      }

      /**
       * Enable all exceptions.
       */
      void enableAllExceptions()
      {
        m_exceptions = ~0;
      }

    private:
      Callback &m_callback; //!< callback for handling events

      std::string m_str; //!< the SMILES/SMARTS string
      std::size_t m_pos; //!< current position in m_str

      // atom
      int m_element; //!< atom element
      int m_isotope; //!< atom isotope
      int m_charge; //!< atom charge
      Chirality m_chiral; //!< atom chirality
      int m_hCount; //!< atom hydrogran count
      int m_class; //!< atom class
      bool m_aromatic; //!< is atom aromatic
      // bond
      int m_bondOrder; //!< last bond's order (1, 2, 3, 4 and 5 for aromatic
      bool m_isUp; //!< is last bond '/'
      bool m_isDown; //!< is last bond '\'
      bool m_explicitBond; //!< is last bond explicit
      // branches & other state
      std::vector<BranchInfo> m_branches; //!< used as stack
      std::map<int, std::vector<RingBondInfo> > m_ringBonds; //!< atom index -> RingBondInfo
      std::vector<ChiralInfo> m_chiralInfo; //!< atom index -> ChiralInfo
      int m_index; //!< current atom index
      int m_prev; // !< previous atom index (-1 for none)
      int m_exceptions;
  };

}

#endif
