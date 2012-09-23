#ifndef SMILEY_MOLECULE_H
#define SMILEY_MOLECULE_H

namespace Smiley {

  struct Atom
  {
    int element;
    int isotope;
    int hCount;
    int charge;
    int atomClass;
    bool aromatic;
  };

  struct Bond
  {
    int source;
    int target;
    int order;
    bool isUp;
    bool isDown;
  };

  struct Molecule
  {
    const Bond* bond(int source, int target) const
    {
      for (std::size_t i = 0; i < bonds.size(); ++i)
        if ((source == bonds[i].source && target == bonds[i].target) || (source == bonds[i].target && target == bonds[i].source))
          return &bonds[i];
      return 0;
    }

    std::vector<Atom> atoms;
    std::vector<Bond> bonds;
    std::map<int, std::pair<Chirality, std::vector<int> > > chirality;
  };

  struct MoleculeSmilesCallback : public CallbackBase
  {
    void clear()
    {
      molecule.atoms.clear();
      molecule.bonds.clear();
      molecule.chirality.clear();
    }

    void addAtom(int atomicNumber, bool aromatic, int isotope, int hCount, int charge, int atomClass)
    {
      molecule.atoms.resize(molecule.atoms.size() + 1);
      molecule.atoms.back().element = atomicNumber;
      molecule.atoms.back().aromatic = aromatic;
      molecule.atoms.back().isotope = isotope;
      molecule.atoms.back().hCount = hCount;
      molecule.atoms.back().charge = charge;
      molecule.atoms.back().atomClass = atomClass;
    }

    void addBond(int source, int target, int order, bool isUp, bool isDown)
    {
      molecule.bonds.resize(molecule.bonds.size() + 1);
      molecule.bonds.back().source = source;
      molecule.bonds.back().target = target;
      molecule.bonds.back().order = order;
      molecule.bonds.back().isUp = isUp;
      molecule.bonds.back().isDown = isDown;
    }

    void setChiral(int index, Chirality chirality, const std::vector<int> &chiralRefs)
    {
      molecule.chirality[index] = std::make_pair(chirality, chiralRefs);
    }

    Molecule molecule;
  };


}

#endif
