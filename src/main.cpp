#include <iostream>
#include <bitset>

using namespace std;

template<int BITS> class uint_kmer {
   using int_type =
      typename std::conditional< BITS <= 8,  uint8_t,
      typename std::conditional< BITS <= 16, uint16_t,
      typename std::conditional< BITS <= 32, uint32_t, uint64_t >::type >::type >::type;
   public:
   int_type i;

   uint_kmer(int a) : i(a){}

   size_t get() {return i;}
};

int main(){
    cout << "hello, i do nothing for now\n";
    
    uint64_t a = 1ul;

    for (int k=1; k<33;){
        a = ((1ul << (2 * (k-1))) - 1);
        cout << k << " " << bitset<64>(a) << " " << a << endl;
        k = k*2;
    }

    cout << sizeof(uint8_t) << " " << sizeof(uint16_t) << " " << sizeof(uint32_t) << " " << sizeof(uint64_t) << endl;
        
    uint_kmer<1> b1 = ((1ul << (20)) - 1); 
    uint_kmer<15> b2 = ((1ul << (20)) - 1);
    uint_kmer<16> b3 = ((1ul << (20)) - 1); 
    uint_kmer<17> b4 = ((1ul << (20)) - 1);
    cout << sizeof(b1) <<" "<< sizeof(b2) <<" "<< sizeof(b3) <<" "<< sizeof(b4) <<" ";
    cout << endl;
    cout << b1.get() <<" "<< b2.get() <<" "<< b4.get();
    cout << endl;
}