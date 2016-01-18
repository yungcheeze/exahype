#include <iostream>

using namespace std;

int main() {
    {
        const int imax=2;
        const int jmax=4;
        double arr1D[imax*jmax];

        for(int i=0;i<imax*jmax;i++) {
            arr1D[i] = i;
        }

        typedef double row_t[jmax];
        row_t *matrix = (row_t *)arr1D;

        for(int i=0;i<imax;i++){
            for(int j=0;j<jmax;j++){
                cout << matrix[i][j] << " ";
            }
        }
        cout << endl << endl;
    }

    //-------

    {
        const int imax=2;
        const int jmax=4;
        const int kmax=3;
        double arr1D[imax*jmax*kmax];

        for(int i=0;i<imax*jmax*kmax;i++) {
            arr1D[i] = i;
        }
        typedef double tensor_t[jmax][kmax]; // 4 3
        tensor_t *tensor = (tensor_t *)arr1D;

        // page - row - column
        cout << tensor[1][3][0] << endl;

        for(int i=0;i<imax;i++){
            for(int j=0;j<jmax;j++){
                for(int k=0;k<kmax;k++){
                    cout << tensor[i][j][k] << " ";
                }
            }
        }

        cout << endl << endl;
    }

    //-----

    {
        const int imax = 2;
        const int jmax = 2;
        const int kmax = 4;
        const int lmax = 4;

        double arr1D[imax*jmax*kmax*lmax];

        for(int i=0;i<imax*jmax*kmax*lmax;i++) {
            arr1D[i] = i;
        }

        typedef double tensor_t[jmax][kmax][lmax];
        tensor_t *tensor = (tensor_t *)arr1D;

        for(int i=0;i<imax;i++){
            for(int j=0;j<jmax;j++){
                for(int k=0;k<kmax;k++){
                    for(int l=0;l<lmax;l++){
                        cout << tensor[i][j][k][l] << " " ;
                    }
                }
            }
        }

    }

    return 0;
}


