#include "smiley.h"
#include "test.h"

using namespace Smiley;

int main(int argc, char **argv)
{
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <smiles string>" << std::endl;
    return 1;  
  }

  PrintCallback callback;
  Parser<PrintCallback> parser(callback);

  try {
    parser.parse(argv[1]);
  } catch (Exception &e) {
    if (e.type() == Exception::SyntaxError)
      std::cerr << "Syntax";
    else
      std::cerr << "Semantics";
    std::cout << "Error: " << e.what() << "." << std::endl;
    std::cout << argv[1] << std::endl;
    for (std::size_t i = 0; i < e.pos(); ++i)
      std::cout << " ";
    for (std::size_t i = 0; i < e.length(); ++i)
      std::cout << "^";
    std::cout << std::endl;
  }
}
