
extern "C" {
  void writeit_();
  void writethat_(int*, double*, char*);
}

int main() {
  writeit_();
  int a(4);
  double b(12.6567676974);
  char c('g');
  writethat_(&a, &b, &c);
}
