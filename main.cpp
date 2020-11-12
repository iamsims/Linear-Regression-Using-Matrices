//minimizing the residual sum of squares for a overdetermined set of equations 
//or linear regression 
//Ax = b: x = ((A`A)^-1)*A*b

#include<bits/stdc++.h> 
using namespace std; 


template <size_t N>
void getCofactor(float (&A)[N][N], float (&temp)[N][N], int p, int q, int n) 
{ 
    int i = 0, j = 0; 
  
    for (int row = 0; row < n; row++) 
    { 
        for (int col = 0; col < n; col++) 
        { 
            if (row != p && col != q) 
            { 
                temp[i][j++] = A[row][col]; 

                if (j == n - 1) 
                { 
                    j = 0; 
                    i++; 
                } 
            } 
        } 
    } 
} 


template <size_t N>
float determinant(float (&A)[N][N], int n) 
{ 
    float D = 0; 

    if (n == 1) 
        return A[0][0]; 
  
    float temp[N][N];  
  
    int sign = 1; 

    for (int f = 0; f < n; f++) 
    { 
        
        getCofactor(A, temp, 0, f, n); 
        D += sign * A[0][f] * determinant(temp, n - 1); 

        sign = -sign; 
    } 
  
    return D; 
} 

template <size_t N>
void adjoint(float (&A)[N][N],float (&adj)[N][N]) 
{ 
    if (N == 1) 
    { 
        adj[0][0] = 1; 
        return; 
    } 
  
    int sign = 1;
    float temp[N][N]; 
  
    for (size_t i=0; i<N; i++) 
    { 
        for (size_t j=0; j<N; j++) 
        { 
         
            getCofactor(A, temp, i, j, N); 
  
            sign = ((i+j)%2==0)? 1: -1; 
  
            adj[j][i] = (sign)*(determinant(temp, N-1)); 
        } 
    } 
}

template <size_t row, size_t col>
void transpose(float (&mat1)[row][col], float (&mat2)[col][row]) {
    for (size_t i = 0; i < row; ++i)
      for (size_t j = 0; j < col; ++j) {
        mat2[j][i] = mat1[i][j];
      }
}

template <size_t N>
void inverse(float (&A)[N][N], float (&inverse)[N][N])
{ 
    // Find determinant of A[][] 
    float det = determinant(A, N);  
 
    float adj[N][N]; 
    adjoint(A, adj); 
  
   
    for (size_t i=0; i<N; i++) 
        for (size_t j=0; j<N; j++) 
            inverse[i][j] = adj[i][j]/float(det); 
  
} 
  

template <size_t row, size_t col>
void display(float (&A)[row][col]) 
{ 
    for (size_t i=0; i<row; i++) 
    { 
        for (size_t j=0; j<col; j++) 
            cout << A[i][j] << " "; 
        cout << endl; 
    } 
} 


template <size_t row1, size_t col1, size_t col2>
void multiply(float (&mat1)[row1][col1], float (&mat2)[col1][col2], float(&res)[row1][col2])
{
    size_t i,j,x;
    for (i = 0; i < row1; i++) 
    {
        for (j = 0; j < col2; j++) 
        {
            res[i][j] = 0;
            for (x = 0; x < col1; x++) 
            {
                *(*(res + i) + j) += *(*(mat1 + i) + x)
                                     * *(*(mat2 + x) + j);
            }
        }
    }
}

// Driver program 
int main() 
{ 
    //HERE: code for inputting the linear equations and parsing

    float A[5][2] = { {0.25, 0.03}, 
                        {-0.37, 0.17},
                        {1.17,1 },
                        {-1.69, -0.17},
                        {1, 1.19}
                    }; 
    float b[5][1]= {{0.97},{0.48},{2.20},{-1.19}, {1.73}};
    float transA[2][5];
    transpose(A,transA);
    //display(transA);
    float AtA[2][2];
    multiply(transA, A, AtA);
    display(AtA);
    float AtAinv[2][2];
    inverse(AtA, AtAinv);
    display(AtAinv);
    float temp[2][5];
    multiply(AtAinv,transA, temp);
    display(temp);
    float ans[2][1];
    multiply(temp, b, ans);
    cout<<"The answer is :"<<endl;
    display(ans);

    return 0; 
} 
