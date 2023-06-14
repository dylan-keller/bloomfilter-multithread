#include <iostream>
#include <bitset>

using namespace std;

int main(){
    cout << "hello, i do nothing for now\n";
    
    uint64_t a = 1ul;

    for (int k=1; k<33; k++){
        a = ((1ul << (2 * (k-1))) - 1);
        cout << k << " " << bitset<64>(a) << " " << a << endl;
    }
}