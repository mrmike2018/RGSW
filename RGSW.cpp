/*
This is an implementation of the RGSW (Ring GSW) scheme using the NTL and GMP C++ libraries.

Sources:
NTL C++ library: https://libntl.org/
GMP C++ library: https://gmplib.org/

Articles and references:
GSW scheme: https://eprint.iacr.org/2013/340.pdf

other related articles
https://web.eecs.umich.edu/~cpeikert/pubs/polyboot.pdf
https://eprint.iacr.org/2021/691.pdf
https://eprint.iacr.org/2020/086.pdf
*/

#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZX.h>
#include <NTL/matrix.h>
#include <NTL/vector.h>


using namespace std;
using namespace NTL;

// public parameters:
// best parameters:
//      bit_length = 32 and base = 256.   or
//      bit_length = 64 and base = 256*256.
ZZ p; // this is the modulous for the ring Z_p.
const long bit_length = 32;
long base = 256; // For now, the implementation supports the base 2 or 256, or 65536.
long K;          // K is defined as bit_length / log2(base)
long poly_degree = 1024; // this is the degree on the modulous polynomial (i.e., n in X^n + 1).

Vec<ZZ_pX> GenerateSecretKey() {
    Vec<ZZ_pX> sk;
    sk.SetLength(2);

    for (size_t i = 0; i < 1; i++) {
        random(sk[i], poly_degree);
    }
    sk[1] = ZZ_pX(1);

    return sk;
}

Mat<ZZ_pX> CreateMatrixA(Vec<ZZ_pX> sk) {

    Mat<ZZ_pX> A;
    A.SetDims(2, 2 * K);

    ZZ_pX x_n_plus_1_poly;
    x_n_plus_1_poly.SetLength(poly_degree + 1);
    x_n_plus_1_poly[0] = 1; x_n_plus_1_poly[poly_degree] = 1;

    // generating a random error vector:
    Vec<ZZ_p> error_vec;
    error_vec.SetLength(2 * K);

    // setting the modulous to 2, for generating a small error vector:
    ZZ_p::init(ZZ(2));

    for (size_t i = 0; i < 2 * K; i++) {
        random(error_vec[i]);
    }

    // setting the modulous back to p:
    ZZ_p::init(p);

    //cout << "\nThe error vector: " << endl;
    //cout << error_vec << endl;

    // in the following lines, we create the matrix A as it is constructed in the GSW scheme:
    Mat<ZZ_pX> A_temp;
    A_temp.SetDims(1, 2 * K);

    ZZ_pX temp_poly;
    for (int i = 0; i < 1; i++) {
        for (int j = 0; j < 2 * K; j++) {
            random(temp_poly, poly_degree);
            A_temp[i][j] = temp_poly;
            A[i][j] = temp_poly;
        }
    }

    // intitializing s_multipliedBy_A_temp to zero
    Vec<ZZ_pX> s_multipliedBy_A_temp;
    s_multipliedBy_A_temp.SetLength(A_temp.NumCols());
    for (size_t i = 0; i < A_temp.NumCols(); i++) {
        s_multipliedBy_A_temp[i] = ZZ_pX(0);
    }

    ZZ_pX FFTmult_result;
    for (size_t i = 0; i < sk.length() - 1; i++) {
        for (size_t j = 0; j < A_temp.NumCols(); j++) {
            for (size_t k = 0; k < A_temp.NumRows(); k++) {
                mul(FFTmult_result, sk[k], A_temp[k][j]);
                s_multipliedBy_A_temp[j] += FFTmult_result;
            }
        }
    }

    for (size_t i = 0; i < s_multipliedBy_A_temp.length(); i++) {
        rem(s_multipliedBy_A_temp[i], s_multipliedBy_A_temp[i], x_n_plus_1_poly);
    }

    // creating matrix A as: e - s * A_temp:
    for (int i = 1; i < 2; i++) {
        for (int j = 0; j < 2 * K; j++) {
            A[i][j] = error_vec[j] - s_multipliedBy_A_temp[j];
        }
    }

    //cout << "matrix A: " << endl;
    //cout << A << endl << endl;

    return A;
}

Mat<ZZ_p> createMatrixG() {

    // creating matrix G:
    Mat<ZZ_p> G;

    G.SetDims(2, 2 * K);
    ZZ_p temp = ZZ_p(1);

    int n = 2; // n: number of rows in the matrix G.
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < 2 * K; j++) {
            G[i][j] = ZZ_p(0);
        }
    }

    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < K; j++) {
            G[i][i * K + j] = temp;
            mul(temp, temp, ZZ_p(base));
        }
        temp = ZZ_p(1);
    }

    //cout << "matrix G: " << endl;
    //cout << G << endl << endl;

    return G;
}

Vec<ZZ> convertZZtoBaseB(ZZ num, long base) {
    unsigned char* pp = new unsigned char[NumBytes(num)];
    BytesFromZZ(pp, num, NumBytes(num));

    Vec<ZZ> result;
    result.SetLength(K);
    long blockSize = log2(base) / 8;
    long blocksToBeRead = blockSize, blocksReadSoFar = 0;

    // for base 256:
    blocksToBeRead = NumBytes(num);
    blocksReadSoFar = 0;
    for (size_t i = 0; i < NumBytes(num); i = i + blockSize) {
        if (NumBytes(num) - blocksReadSoFar >= blockSize) {
            blocksToBeRead = blockSize;
        }
        else {
            blocksToBeRead = NumBytes(num) - blocksReadSoFar;
        }

        ZZFromBytes(result[i / blockSize], pp + i, blocksToBeRead);
        blocksReadSoFar = blocksReadSoFar + blockSize;
    }

    return result;
}

Vec<ZZX> G_inverse(ZZ_pX polynomial) {
    // this function gets a polynomial of degree at most 'poly_degree - 1', in which
    // each coefficient is a number of length 'bit_length'.
    // Then the function decompses each coefficients into its binary representation, 
    // and forms a vector of polynomials with binary coefficients the way that G^-1 in the GSW scheme works:

    // variables for converting to an arbitrary base:
    ZZX poly_with_baseB_coeffs;
    Vec<ZZX> coeffs_convertedToBaseB;
    poly_with_baseB_coeffs.SetLength(poly_degree);
    coeffs_convertedToBaseB.SetLength(K);

    long padding_length = poly_degree - polynomial.rep.length();
    if (polynomial.rep.length() < poly_degree) {
        for (size_t i = 0; i < padding_length; i++) {
            polynomial.rep.append(ZZ_p(0));
        }
    }

    ZZ temp_coeff;
    Vec<ZZ> conversionToBaseBResult;
    
    // we convert the input polynomial, which is of type ZZ_pX, to a poly of type ZZX
    // so that we can easily convert its coefficients to baseB
    ZZX polynomial_as_ZZX;
    polynomial_as_ZZX.SetLength(poly_degree);
    conv(polynomial_as_ZZX, polynomial);

    // padding the polynomial:
    long padding_length2 = poly_degree - polynomial_as_ZZX.rep.length();
    if (polynomial_as_ZZX.rep.length() < poly_degree) {
        for (size_t i = 0; i < padding_length2; i++) {
            polynomial_as_ZZX.rep.append(ZZ(0));
        }
    }

    // this is only for handling base 2:
    if (base == 2) {
        for (size_t i = 0; i < K; i++) {
            for (size_t j = 0; j < poly_degree; j++) {
                poly_with_baseB_coeffs[j] = bit(polynomial_as_ZZX[j], i);
            }
            coeffs_convertedToBaseB[i] = poly_with_baseB_coeffs;
        }

        return coeffs_convertedToBaseB;
    }

    // initializing all elements of coeffs_convertedToBaseB to zero
    for (size_t i = 0; i < K; i++) {
        coeffs_convertedToBaseB[i] = ZZX(0);
    }

    // this is for base 256 or other bases:
    Vec<ZZ> coeff_j_convertedToBaseB;
    for (size_t j = 0; j < poly_degree; j++) {
        coeff_j_convertedToBaseB = convertZZtoBaseB(polynomial_as_ZZX[j], base);
        for (size_t i = 0; i < K; i++) {
            coeffs_convertedToBaseB[i].rep.append(coeff_j_convertedToBaseB[i]);
        }
    }

    return coeffs_convertedToBaseB;
}

Mat<ZZX> G_inverse_of_ciphertext(Mat<ZZ_pX> ciphertext) {
    // this function  gets a ciphertext as input, which is a matrix of dimension (2 * 2*bit_length).
    // The function returns the bit-decomposition or byte-decomposition of the ciphertext.
    // the output is a matrix of dimension (2*bit_length * 2*bit_length).
    // each element of the result matrix is a polynomial with binary coefficients (when base is 2),
    // or a polynomial with coefficients {0, 1, 2, ..., 7} when the base is set to 256.

    Mat<ZZX> ciphertext_decomposedInBaseB;
    ciphertext_decomposedInBaseB.SetDims(2 * K, 2 * K);

    Vec<ZZX> poly_coeffs_decomposed_inBaseB;
    poly_coeffs_decomposed_inBaseB.SetLength(K);

    for (size_t i = 0; i < 2 * K; i++) {
        for (size_t j = 0; j < 2; j++) {
            poly_coeffs_decomposed_inBaseB = G_inverse(ciphertext[j][i]);
            for (size_t k = 0; k < K; k++) {
                ciphertext_decomposedInBaseB[j * K + k][i] = poly_coeffs_decomposed_inBaseB[k];
            }
        }
    }

    //cout << "\nciphertext_decomposedInBaseB: " << endl;
    //cout << ciphertext_decomposedInBaseB << endl;

    return ciphertext_decomposedInBaseB;
}

// Implementation of the encryption algorithm for the RGSW scheme
Mat<ZZ_pX> encrypt(int message, Vec<ZZ_pX> sk) {

    ZZ_pX x_n_plus_1_poly;
    x_n_plus_1_poly.SetLength(poly_degree + 1);
    x_n_plus_1_poly[0] = 1; x_n_plus_1_poly[poly_degree] = 1;

    ZZ_p m;
    m = ZZ_p(message);

    // creating matrix A:
    Mat<ZZ_pX> A;
    A = CreateMatrixA(sk);

    //cout << "matrix A: " << endl;
    //cout << A << endl << endl;

    // constructing the random matrix R with binary entries (as defined in the GSW scheme):
    ZZ_p::init(ZZ(2)); // setting the modulous to two for creating random polynomials with small coefficients.
    Mat<ZZ_pX> R;
    R.SetDims(2 * K, 2 * K);

    for (int i = 0; i < 2 * K; i++) {
        for (int j = 0; j < 2 * K; j++) {
            random(R[i][j], poly_degree);
        }
    }

    //cout << "matrix R: " << endl;
    //cout << R << endl << endl;

    // setting the modulous back to p:
    ZZ_p::init(p);

    // calculating A * R:
    Mat<ZZ_pX> A_milttipliedBy_R;
    A_milttipliedBy_R.SetDims(2, 2 * K);

    ZZ_pX FFTmult_result;
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2 * K; j++) {
            for (int k = 0; k < 2 * K; k++) {
                mul(FFTmult_result, A[i][k], R[k][j]);
                add(A_milttipliedBy_R[i][j], A_milttipliedBy_R[i][j], FFTmult_result);
            }
        }
    }

    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2 * K; j++) {
            rem(A_milttipliedBy_R[i][j], A_milttipliedBy_R[i][j], x_n_plus_1_poly);
        }
    }

    // creating matrix G:
    Mat<ZZ_p> G;
    G = createMatrixG();

    //cout << "matrix G: " << endl;
    //cout << G << endl << endl;

    // calculating m * G:
    Mat<ZZ_p> m_milttipliedBy_G;
    m_milttipliedBy_G.SetDims(2, 2 * K);

    if (message == 1)
        m_milttipliedBy_G = G;
    else {
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2 * K; j++) {
                m_milttipliedBy_G[i][j] = ZZ_p(0);
            }
        }
    }

    //cout << "m_milttipliedBy_G: " << endl;
    //cout << m_milttipliedBy_G << endl;

    // Encrypting a message (where m can be 0 or 1) using the RGSW scheme:
    // ciphertext or Enc(m) = A . R + m . G:
    Mat<ZZ_pX> ciphertext;
    ciphertext.SetDims(2, 2 * K);

    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2 * K; j++) {
            ciphertext[i][j] = A_milttipliedBy_R[i][j] + m_milttipliedBy_G[i][j];
        }
    }

    //cout << "\nciphertext: \n" << ciphertext << endl;

    return ciphertext;
}

ZZ_pX dot_product(Vec<ZZ_pX> v1, Vec<ZZ_pX> v2) {

    ZZ_pX x_n_plus_1_poly;
    x_n_plus_1_poly.SetLength(poly_degree + 1);
    x_n_plus_1_poly[0] = 1; x_n_plus_1_poly[poly_degree] = 1;

    if (v1.length() != v2.length())
        cout << "The lengths of vectors are different!" << endl;

    ZZ_pX result = ZZ_pX(0);
    ZZ_pX FFTmult_result;
    for (size_t i = 0; i < v1.length(); i++) {
        mul(FFTmult_result, v1[i], v2[i]);
        result += FFTmult_result;
    }

    rem(result, result, x_n_plus_1_poly);

    return result;
}

Vec<ZZ_pX> product_vector_matrix(Vec<ZZ_pX> vec, Mat<ZZ_pX> mat) {

    ZZ_pX x_n_plus_1_poly;
    x_n_plus_1_poly.SetLength(poly_degree + 1);
    x_n_plus_1_poly[0] = 1; x_n_plus_1_poly[poly_degree] = 1;

    if (vec.length() != mat.NumRows()) {
        cout << "Error: the dimnesions of the vector and matrix mismatch!" << endl;
        exit(-1);
    }

    Vec<ZZ_pX> result;
    result.SetLength(mat.NumCols());

    for (size_t i = 0; i < mat.NumCols(); i++) {
        result[i] = ZZ_pX(0);
    }

    ZZ_pX FFTmult_result;
    for (size_t i = 0; i < mat.NumCols(); i++) {
        for (size_t j = 0; j < vec.length(); j++) {
            mul(FFTmult_result, vec[j], mat[j][i]);
            result[i] += FFTmult_result;
        }
    }

    for (size_t i = 0; i < result.length(); i++) {
        rem(result[i], result[i], x_n_plus_1_poly);
    }

    return result;
}

// Implementation of the decryption algorithm for the RGSW scheme:
// The following decryption algorithm works better for base = 2.
// The algorithm checks whether the result falls in (-p/4, p/4) or in (p/4, p/2) or (-p/2, -p/4)
//int decrypt(Vec<ZZ_pX> sk, Mat<ZZ_pX> ciphertext) {
int decrypt_disable_temporarily(Vec<ZZ_pX> sk, Mat<ZZ_pX> ciphertext) {
    int decrypted_message = 0;
    int decrypted_message_with_last_col = 0;

    Vec<ZZ_pX> penultimate_column_of_c, last_column_of_c;
    Vec<ZZ_pX> first_column_of_c, second_column_of_c;

    penultimate_column_of_c.append(ciphertext[0][2 * K - 2]);
    penultimate_column_of_c.append(ciphertext[1][2 * K - 2]);

    last_column_of_c.append(ciphertext[0][2 * K - 1]);
    last_column_of_c.append(ciphertext[1][2 * K - 1]);

    first_column_of_c.append(ciphertext[0][0]);
    first_column_of_c.append(ciphertext[1][0]);

    second_column_of_c.append(ciphertext[0][1]);
    second_column_of_c.append(ciphertext[1][1]);

    ZZ_pX result1, result2, result3, result4;

    result3 = dot_product(sk, penultimate_column_of_c);
    result4 = dot_product(sk, last_column_of_c);

    ZZ result3ZZ, result4ZZ, p_dividedBy_2, p_dividedBy_4;

    // getting the constant term of the polynomial and
    // converting it to type ZZ (type for big integer) to be fed into the NTL's compare() function:
    conv(result3ZZ, result3[0]);
    conv(result4ZZ, result4[0]);
    div(p_dividedBy_2, p, ZZ(2));
    div(p_dividedBy_4, p, ZZ(4));

    // if 0 < result3 < p/4, then decrypt the message to 0.
    // if p/4 < result3 < p/2, then decrypt the message to 1.
    cout << "Note that result3 = dot_product(sk, penultimate_column_of_c)." << endl;

    // if result3ZZ is in (0, p/4), then decrypt to 0
    if (compare(result3ZZ, ZZ(0)) > 0 && compare(result3ZZ, p_dividedBy_4) < 0) {
        // i.e., result3ZZ is in (0, p/4)
        decrypted_message = 0;
    }

    // if result3ZZ is in (p/4, p/2), then decrypt to 1
    if (compare(result3ZZ, p_dividedBy_4) > 0 && compare(result3ZZ, p_dividedBy_2) < 0) {
        // i.e., result3ZZ is in (p/4, p/2)
        decrypted_message = 1;
    }

    // if result3ZZ > p/2, then we have to map it to (-p/2, 0):
    if (compare(result3ZZ, p_dividedBy_2) > 0) {
        ZZ result3ZZ_minus_p = result3ZZ - p;

        // if -p/2 < result3ZZ_minus_p < -p/4, then decrypt to 1:
        if (compare(result3ZZ_minus_p, -p_dividedBy_2) > 0 && compare(result3ZZ_minus_p, -p_dividedBy_4) < 0) {
            decrypted_message = 1;
        }

        // if -p/4 < result3ZZ_minus_p < 0, then again decrypt to 0:
        if (compare(result3ZZ_minus_p, -p_dividedBy_4) > 0 && compare(result3ZZ_minus_p, ZZ(0)) < 0) {
            decrypted_message = 0;
        }
    }

    return decrypted_message;
}

// this decryption algorithm works better with different bases, including base 256
int decrypt(Vec<ZZ_pX> sk, Mat<ZZ_pX> ciphertext) {
    int decrypted_message = 0;

    Vec<ZZ_pX> penultimate_column_of_c, last_column_of_c;

    penultimate_column_of_c.append(ciphertext[0][2 * K - 2]);
    penultimate_column_of_c.append(ciphertext[1][2 * K - 2]);

    last_column_of_c.append(ciphertext[0][2 * K - 1]);
    last_column_of_c.append(ciphertext[1][2 * K - 1]);

    ZZ_pX result1, result2, result3, result4;
    result3 = dot_product(sk, penultimate_column_of_c);
    result4 = dot_product(sk, last_column_of_c);

    // this means to use that last column for decryption:
    // it seems that when the base is 256, using the last column works better for decryption.
    if (base == 256) {
        result3 = result4;
    }

    ZZ_p constant_term_of_result3, constant_term_of_result4;

    // getting the constant term of the result3 polynomial:
    if (result3 == ZZ_pX(0)) {
        constant_term_of_result3 = ZZ_p(0);
    }
    else {
        constant_term_of_result3 = result3[0];
    }
    // getting the constant term of the result4 polynomial:
    if (result4 == ZZ_pX(0)) {
        constant_term_of_result4 = ZZ_p(0);
    }
    else {
        constant_term_of_result4 = result4[0];
    }

    ZZ result3ZZ, result4ZZ;
    long comparison_result3, comparison_result4;

    // converting constant_term_of_result3 to type ZZ to be fed into the NTL's compare() function:
    conv(result3ZZ, constant_term_of_result3);

    ZZ p_dividedBy_2, p_dividedBy_4;
    div(p_dividedBy_2, p, ZZ(2));
    div(p_dividedBy_4, p, ZZ(4));

    ZZ base_to_power_K_minus_two, base_to_power_K_minus_one, base_to_power_K;
    if (K < 2) {
        cout << "Error: The selected parameteres and/or dimensions of the RGSW are too small!" << endl;
        exit(-1);
    }
    power(base_to_power_K_minus_two, ZZ(base), K - 2);
    power(base_to_power_K_minus_one, ZZ(base), K - 1);
    power(base_to_power_K, ZZ(base), K);

    ZZ two_to_power_l_minus_two, two_to_power_l_minus_one;
    power(two_to_power_l_minus_two, ZZ(2), bit_length - 2);
    power(two_to_power_l_minus_one, ZZ(2), bit_length - 1);

    ZZ distance_to_zero, distance_to_prime_p;
    ZZ distance_to_base_to_K_minus_two;
    ZZ distance_to_two_to_l_minus_one, distance_to_two_to_l_minus_two;

    distance_to_zero = abs(result3ZZ - ZZ(0));
    distance_to_prime_p = abs(result3ZZ - p);
    distance_to_two_to_l_minus_two = abs(result3ZZ - two_to_power_l_minus_two);
    distance_to_two_to_l_minus_one = abs(result3ZZ - two_to_power_l_minus_one);

    distance_to_base_to_K_minus_two = abs(result3ZZ - base_to_power_K_minus_two);

    // in this decryption approach we basically check if result3ZZ
    // is closer to zero or prime or
    // closer to base_to_K_minus_two.
    // in the first case, we decrypt to zero, in the second case we decrypt to one.

    if ((result3ZZ < base_to_power_K_minus_two) &&
        (distance_to_zero < distance_to_base_to_K_minus_two)) {
        decrypted_message = 0;
    }

    if ((result3ZZ < base_to_power_K_minus_two) &&
        (distance_to_base_to_K_minus_two < distance_to_zero)) {
        decrypted_message = 1;
    }

    if ((result3ZZ >= base_to_power_K_minus_two) &&
        (distance_to_base_to_K_minus_two < distance_to_prime_p)) {
        decrypted_message = 1;
    }

    if ((result3ZZ > base_to_power_K_minus_two) &&
        (distance_to_prime_p < distance_to_base_to_K_minus_two)) {
        decrypted_message = 0;
    }

    return decrypted_message;
}

// Implementation of the addition gate for the RGSW scheme
Mat<ZZ_pX> add(Mat<ZZ_pX> c1, Mat<ZZ_pX> c2) {
    Mat<ZZ_pX> addtion_result;

    if ((c1.NumRows() != c2.NumRows()) || (c1.NumCols() != c2.NumCols())) {
        cout << "Error in the 'add' function: the dimensions of the ciphertext do not matach!" << endl;
        exit(-1);
    }

    addtion_result.SetDims(c1.NumRows(), c1.NumCols());

    for (size_t i = 0; i < c1.NumRows(); i++) {
        for (size_t j = 0; j < c1.NumCols(); j++) {
            addtion_result[i][j] = c1[i][j] + c2[i][j];
        }
    }

    //cout << "\naddtion_result: \n" << addtion_result << endl;

    return addtion_result;
}

// Implementation of the multiplication gate for the RGSW scheme:
Mat<ZZ_pX> multiply(Mat<ZZ_pX> c1, Mat<ZZ_pX> c2) {

    ZZ_pX x_n_plus_1_poly;
    x_n_plus_1_poly.SetLength(poly_degree + 1);
    x_n_plus_1_poly[0] = 1; x_n_plus_1_poly[poly_degree] = 1;

    if ((c1.NumRows() != c2.NumRows()) || (c1.NumCols() != c2.NumCols())) {
        cout << "multiply function: the dimensions of the ciphertext do not matach!" << endl;
        exit(-1);
    }

    Mat<ZZ_pX> multiplication_result;
    multiplication_result.SetDims(c1.NumRows(), c2.NumCols());

    Mat<ZZX> G_inverseOf_c2;
    G_inverseOf_c2.SetDims(2 * K, 2 * K);

    // calculating G^-1(c2):
    G_inverseOf_c2 = G_inverse_of_ciphertext(c2);

    //cout << "G_inverseOf_c2: \n" << G_inverseOf_c2 << endl << endl;    

    Mat<ZZ_pX> G_inverseOf_c2_ZZ_pX;
    G_inverseOf_c2_ZZ_pX.SetDims(G_inverseOf_c2.NumRows(), G_inverseOf_c2.NumCols());

    // converting G_inverse_of_ciphertext from ZZX to ZZ_pX:
    for (size_t i = 0; i < G_inverseOf_c2.NumRows(); i++) {
        for (size_t j = 0; j < G_inverseOf_c2.NumCols(); j++) {
            conv(G_inverseOf_c2_ZZ_pX[i][j], G_inverseOf_c2[i][j]);
        }
    }

    // initializing all elements of multiplication_result to zero.
    multiplication_result.SetDims(c1.NumRows(), c2.NumCols());
    for (size_t i = 0; i < c1.NumRows(); i++) {
        for (size_t j = 0; j < c2.NumCols(); j++) {
            multiplication_result[i][j] = ZZ_pX(0);
        }
    }

    ZZ_pX FFTmult_result;
    for (size_t i = 0; i < c1.NumRows(); i++) {
        for (size_t j = 0; j < c1.NumCols(); j++) {
            for (size_t k = 0; k < G_inverseOf_c2_ZZ_pX.NumRows(); k++) {
                mul(FFTmult_result, c1[i][k], G_inverseOf_c2_ZZ_pX[k][j]);
                add(multiplication_result[i][j], multiplication_result[i][j], FFTmult_result);
            }
        }
    }

    for (size_t i = 0; i < c1.NumRows(); i++) {
        for (size_t j = 0; j < c1.NumCols(); j++) {
            rem(multiplication_result[i][j], multiplication_result[i][j], x_n_plus_1_poly);
        }
    }

    //cout << "\nmultiplication_result: \n" << multiplication_result << endl;

    return multiplication_result;
}

void RGSW() {

    GenPrime(p, bit_length);
    cout << "The modulus (i.e., prime p) is: " << p;
    cout << "   (it is a " << bit_length << " bits prime number)" << endl;
    cout << "The base is: " << base << endl;
    ZZ_p::init(p);

    if (base != 2 && base != 256 && base != 65536) {
        cout << "Fow now, homomorphic multiplication works only with base 2 or 256." << endl;
        cout << "Please change the base to either 2 or 256 or disable the the multiply(c1, c2) function." << endl;
        cout << "and run the program again!" << endl;

        exit(0);
    }

    K = bit_length / log2(base);
    if (K < 1) {
        K = 1;
    }

    cout << "K := bit_length / log2(base) is: " << K << endl;

    cout << "The degree of the modulus polynomial (i.e., n in x^n + 1) is: " << poly_degree << endl << endl;
    cout << "Please wait......., some computations might take a while ..." << endl;

    ///////////////////////////////// generating the secret key:
    Vec<ZZ_pX> sk;
    sk = GenerateSecretKey();

    ///////////////////////////////// initializing messages and ciphertexts:
    int message1 = 1, message2 = 1, message3 = 1, message4 = 0;
    Mat<ZZ_pX> ciphertext1, ciphertext2, addition_result_ct, multiplication_result_ct;
    int decryption_result1, decryption_result2;
    int decryption_result_addition, decryption_result_multiplication;
    double t; // for measurig time

    cout << "Message1 = " << message1 << ", and Message2 = " << message2 << endl;

    ciphertext1 = encrypt(message1, sk);

    ciphertext2 = encrypt(message2, sk);

    ///////////////////////////////// decrypting the result:
    decryption_result1 = decrypt(sk, ciphertext1);
    decryption_result2 = decrypt(sk, ciphertext2);

    cout << "\nCiphertext1 was decrypted to " << decryption_result1 << ". ";
    cout << "\nCiphertext2 was decrypted to " << decryption_result2 << ". ";

    ///////////////////////////////// Performing homomorphic operations:
    addition_result_ct = add(ciphertext1, ciphertext2);

    t = GetTime();
    multiplication_result_ct = multiply(ciphertext1, ciphertext2);
    t = GetTime() - t;
    cout << "\n\nMultiplying two RGSW ciphertexts took: ";
    cout << t << " (seconds)" << endl;

    decryption_result_addition = decrypt(sk, addition_result_ct);
    decryption_result_multiplication = decrypt(sk, multiplication_result_ct);

    cout << "\nDecryption of the addition of the ciphertexts: "
        << decryption_result_addition << endl;    

    cout << "Decryption of the multiplication of the ciphertexts: "
        << decryption_result_multiplication << endl;

}


int main() {

    cout << "################################################################" << endl;
    cout << "#------------- Welcome to RGSW Implementation! -------------#" << endl;
    cout << "################################################################" << endl;

    RGSW();

    cout << endl;
    cout << "################################################################" << endl;
    cout << "#----------------------- End of Program -----------------------#" << endl;
    cout << "################################################################" << endl;

    return 0;
}