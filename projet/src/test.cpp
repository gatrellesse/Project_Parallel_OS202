#include <bits/stdc++.h>

using namespace std;

int main(){

    vector<pair<int, int>> v = {{1, 2}, {3, 4}, {5, 6}};

    for (auto &[a, b] : v){
        cout << a << " " << b << endl;
        b = 0;
    }
    for (auto &[a, b] : v){
        cout << a << " " << b << endl;
        // b = 0;
    }
}