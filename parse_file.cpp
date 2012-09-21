#include "smiley.h"
#include "molecule.h"

#include <fstream>

using namespace Smiley;

int main(int argc, char **argv)
{
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <smiles file>" << std::endl;
    return 1;  
  }

  MoleculeSmilesCallback callback;
  Parser<MoleculeSmilesCallback> parser(callback);

  std::ifstream ifs(argv[1]);
  std::string line;
  while (std::getline(ifs, line)) {
    callback.molecule.atoms.clear();
    callback.molecule.bonds.clear();

    try {
      parser.parse(line);
    } catch (Exception &e) {
      if (e.type() == Exception::SyntaxError)
        std::cerr << "Syntax";
      else
        std::cerr << "Semantics";
      std::cerr << "Error: " << e.what() << "." << std::endl;
      std::cerr << line << std::endl;
      for (std::size_t i = 0; i < e.pos(); ++i)
        std::cerr << " ";
      for (std::size_t i = 0; i < e.length(); ++i)
        std::cerr << "^";
      std::cerr << std::endl;
    }
  }
}
