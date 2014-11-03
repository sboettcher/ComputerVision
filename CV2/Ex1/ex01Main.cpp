// Copyright 2014 Sebastian Boettcher

#include <NMath.cpp>
#include <assert.h>
#include <math.h>
#include <vector>
#include "Point2D.cpp"


bool readData(const char* file, std::vector<Point2D>& data) {
  std::ifstream myReadFile;
  myReadFile.open(file);
  float X, Y;
  if (myReadFile.is_open()) {
    while (myReadFile >> X >> Y) {
      data.push_back(Point2D(X, Y));
    }
  } else {
    return false;
  }
  myReadFile.close();
  return true;
}

void fillLinSystemMat(CMatrix<float>& mat, const std::vector<Point2D>& model, const std::vector<Point2D>& data) {
  assert(mat.ySize() == static_cast<int>(2*model.size()));
  assert(mat.ySize() == static_cast<int>(2*data.size()));
  // std::cout << H.xSize() << " " << model.size() << " " << data.size() << std::endl;
  for (int i = 0; i < static_cast<int>(model.size()); ++i) {
    mat(0, 2*i+1) = model[i].getX();
    mat(1, 2*i+1) = model[i].getY();
    mat(2, 2*i+1) = 1.0;
    mat(3, 2*i) = -model[i].getX();
    mat(4, 2*i) = -model[i].getY();
    mat(5, 2*i) = -1.0;
    mat(6, 2*i) = data[i].getY() * model[i].getX();
    mat(6, 2*i+1) = -data[i].getX() * model[i].getX();
    mat(7, 2*i) = data[i].getY() * model[i].getY();
    mat(7, 2*i+1) = -data[i].getX() * model[i].getY();
    mat(8, 2*i) = data[i].getY();
    mat(8, 2*i+1) = -data[i].getX();
  }
}

void fillConstraintMat(CMatrix<float>& mat, const std::vector<CMatrix<float> >& H) {
  assert(mat.xSize() == 6);
  assert(mat.ySize() == static_cast<int>(2*H.size()));
  for (int i = 0; i < static_cast<int>(H.size()); ++i) {
    mat(0, 2*i) = H[i](0, 0) * H[i](1, 0);
    mat(1, 2*i) = H[i](0, 0) * H[i](1, 1) + H[i](0, 1) * H[i](1, 0);
    mat(2, 2*i) = H[i](0, 1) * H[i](1, 1);
    mat(3, 2*i) = H[i](0, 2) * H[i](1, 0) + H[i](0, 0) * H[i](1, 2);
    mat(4, 2*i) = H[i](0, 2) * H[i](1, 1) + H[i](0, 1) * H[i](1, 2);
    mat(5, 2*i) = H[i](0, 2) * H[i](1, 2);

    mat(0, 2*i+1) = H[i](0, 0) * H[i](0, 0) - H[i](1, 0) * H[i](1, 0);
    mat(1, 2*i+1) = H[i](0, 0) * H[i](0, 1) + H[i](0, 1) * H[i](0, 0) - H[i](1, 0) * H[i](1, 1) + H[i](1, 1) * H[i](1, 0);
    mat(2, 2*i+1) = H[i](0, 1) * H[i](0, 1) - H[i](1, 1) * H[i](1, 1);
    mat(3, 2*i+1) = H[i](0, 2) * H[i](0, 0) + H[i](0, 0) * H[i](0, 2) - H[i](1, 2) * H[i](1, 0) + H[i](1, 0) * H[i](1, 2);
    mat(4, 2*i+1) = H[i](0, 2) * H[i](0, 1) + H[i](0, 1) * H[i](0, 2) - H[i](1, 2) * H[i](1, 1) + H[i](1, 1) * H[i](1, 2);
    mat(5, 2*i+1) = H[i](0, 2) * H[i](0, 2) - H[i](1, 2) * H[i](1, 2);
  }
}

void getHomographyMat(const CMatrix<float>& mat, CMatrix<float>& H) {
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      H(j, i) = mat(8, 3*i+j);
    }
  }
}

float getK(CMatrix<float>& K, const CMatrix<float>& B) {
  K(2, 1) = (B(0, 1) * B(0, 2) - B(0, 0) * B(1, 2)) / (B(0, 0) * B(1, 1) - B(0, 1) * B(0, 1));
  float lambda = B(2, 2) - (B(0, 2) * B(0, 2) + K(2, 1)*(B(0, 1) * B(0, 2) - B(0, 0) * B(1, 2))) / B(0, 0);
  K(0, 0) = sqrt(lambda / B(0, 0));
  K(1, 1) = sqrt(lambda * B(0, 0) / (B(0, 0) * B(1, 1) - B(0, 1) * B(0, 1)));
  K(1, 0) = -B(0, 1) * K(0, 0) * K(0, 0) * K(1, 1) / lambda;
  K(2, 0) = K(1, 0) * K(2, 1) / K(1, 1) - B(0, 2) * K(0, 0) * K(0, 0) / lambda;
  return lambda;
}

float getVecNorm(const CVector<float>& vec) {
  float norm = 0.0;
  for (int i = 0; i < vec.size(); ++i) {
    norm += pow(vec(i), 2);
  }
  return sqrt(norm);
}

float getMatNorm(const CMatrix<float>& mat) {
  float norm = 0.0;
  for (int i = 0; i < mat.xSize(); ++i) {
    for (int j = 0; j < mat.ySize(); ++j) {
      norm += pow(mat(i, j), 2);
    }
  }
  return sqrt(norm);
}

void printCVec(const CVector<float>& vec) {
  std::cout << "[";
  for (int i = 0; i < vec.size()-1; ++i) {
    std::cout << vec[i] << ", ";
  }
  // std::cout << vec[vec.size()-1] << "] Norm: " << getVecNorm(vec) << std::endl;
  std::cout << vec[vec.size()-1] << "]" << std::endl;
}

void printCMat(const CMatrix<float>& mat) {
  CVector<float> row(mat.xSize(), 0.0);
  std::cout << "[ Matrix: " << mat.ySize() << "x" << mat.xSize() << std::endl;
  for (int i = 0; i < mat.ySize(); ++i) {
    mat.getVector(row, i);
    printCVec(row);
  }
  std::cout << "] Norm: " << getMatNorm(mat) << std::endl;
}

void matAbs(CMatrix<float>& mat) {
  for (int i = 0; i < mat.xSize(); ++i) {
    for (int j = 0; j < mat.ySize(); ++j) {
      mat(i, j) = fabs(mat(i, j));
    }
  }
}

void normalizePoints(std::vector<Point2D>& data) {
  float sumX = 0.0;
  float sumY = 0.0;
  for (size_t i = 0; i < data.size(); ++i) {
    sumX += data[i].getX();
    sumY += data[i].getY();
  }
  sumX = 1.0/sumX;
  sumY = 1.0/sumY;
  for (size_t i = 0; i < data.size(); ++i) {
    data[i].setX(data[i].getX() * sumX);
    data[i].setY(data[i].getY() * sumY);
  }
}

int main(int argc, const char* argv[]) {
  std::vector<Point2D> model;
  std::vector<Point2D> data1;
  std::vector<Point2D> data2;
  std::vector<Point2D> data3;
  std::vector<Point2D> data4;
  std::vector<Point2D> data5;
  readData("model.txt", model);
  readData("data1.txt", data1);
  readData("data2.txt", data2);
  readData("data3.txt", data3);
  readData("data4.txt", data4);
  readData("data5.txt", data5);

  normalizePoints(model);
  normalizePoints(data1);
  normalizePoints(data2);
  normalizePoints(data3);
  normalizePoints(data4);
  normalizePoints(data5);

  int size = 2 * model.size();
  CMatrix<float> A1(9, size, 0.0);
  CMatrix<float> A2(9, size, 0.0);
  CMatrix<float> A3(9, size, 0.0);
  CMatrix<float> A4(9, size, 0.0);
  CMatrix<float> A5(9, size, 0.0);
  fillLinSystemMat(A1, model, data1);
  fillLinSystemMat(A2, model, data2);
  fillLinSystemMat(A3, model, data3);
  fillLinSystemMat(A4, model, data4);
  fillLinSystemMat(A5, model, data5);

  CMatrix<float> S1(A1.xSize(), A1.xSize(), 0.0);
  CMatrix<float> S2(A2.xSize(), A2.xSize(), 0.0);
  CMatrix<float> S3(A3.xSize(), A3.xSize(), 0.0);
  CMatrix<float> S4(A4.xSize(), A4.xSize(), 0.0);
  CMatrix<float> S5(A5.xSize(), A5.xSize(), 0.0);
  CMatrix<float> V1(A1.xSize(), A1.xSize(), 0.0);
  CMatrix<float> V2(A2.xSize(), A2.xSize(), 0.0);
  CMatrix<float> V3(A3.xSize(), A3.xSize(), 0.0);
  CMatrix<float> V4(A4.xSize(), A4.xSize(), 0.0);
  CMatrix<float> V5(A5.xSize(), A5.xSize(), 0.0);

  NMath::svd(A1, S1, V1, true, 40);
  NMath::svd(A2, S2, V2, true, 40);
  NMath::svd(A3, S3, V3, true, 40);
  NMath::svd(A4, S4, V4, true, 40);
  NMath::svd(A5, S5, V5, true, 40);

  CMatrix<float> H1(3, 3, 0.0);
  CMatrix<float> H2(3, 3, 0.0);
  CMatrix<float> H3(3, 3, 0.0);
  CMatrix<float> H4(3, 3, 0.0);
  CMatrix<float> H5(3, 3, 0.0);
  getHomographyMat(V1, H1);
  getHomographyMat(V2, H2);
  getHomographyMat(V3, H3);
  getHomographyMat(V4, H4);
  getHomographyMat(V5, H5);

  // printCMat(H1);
  // printCMat(H2);
  // printCMat(H3);
  // printCMat(H4);
  // printCMat(H5);

  std::vector<CMatrix<float> > H;
  H.push_back(H1);
  H.push_back(H2);
  H.push_back(H3);
  H.push_back(H4);
  H.push_back(H5);

  CMatrix<float> A(6, 10, 0.0);
  CMatrix<float> S(6, 6, 0.0);
  CMatrix<float> V(6, 6, 0.0);
  CMatrix<float> B(3, 3, 0.0);

  fillConstraintMat(A, H);
  NMath::svd(A, S, V, true, 40);

  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 2; ++j) {
      B(j, i) = V(5, 2*i+j);
    }
  }

  printCMat(B);

  CMatrix<float> K(3, 3, 0.0);
  float lambda = getK(K, B);

  std::cout << "K:" << std::endl;
  printCMat(K);
  std::cout << "lambda: " << lambda << std::endl;
}












