#include <bitset>
#include <iostream>

using namespace std;

int main() {
  for (unsigned int i = 0; i < 20; ++i) {
    bitset<4> foo(i);
    bitset<8> bar(i);
    cout << "" << i << " -> " << foo << ", " << bar << endl;
  }
  return EXIT_SUCCESS;
}
