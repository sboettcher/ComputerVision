// Copyright 2014 Sebastian Boettcher

#include <CMatrix.h>
#include <assert.h>
#include <algorithm>

// draw borders to the edge of the image, with the given pixel strength and color
void drawBorder(CMatrix<float> & image, const int & strength = 1, const float & color = 0) {
  for (int i = 0; i < strength; ++i) {
    image.drawLine(0, i, image.xSize(), i, color);
    image.drawLine(i, 0, i, image.ySize(), color);
    image.drawLine(image.xSize() - i, 0, image.xSize() - i, image.ySize(), color);
    image.drawLine(0, image.ySize() - i, image.xSize(), image.ySize() - i, color);
  }
}

// do Haar decomposition
void decompHaar(const CMatrix<float> & input, CMatrix<float> & output, const int & lvl = 2, const bool & borders = false) {
  assert(input.xSize() == output.xSize());
  assert(input.ySize() == output.ySize());

  CMatrix<float> c_tmp(input);
  CMatrix<float> c;
  CMatrix<float> d_h;
  CMatrix<float> d_v;
  CMatrix<float> d_d;

  for (int l = 1; l <= lvl; ++l) {
    c.setSize(c_tmp.xSize() / 2, c_tmp.ySize() / 2);
    d_h.setSize(c_tmp.xSize() / 2, c_tmp.ySize() / 2);
    d_v.setSize(c_tmp.xSize() / 2, c_tmp.ySize() / 2);
    d_d.setSize(c_tmp.xSize() / 2, c_tmp.ySize() / 2);

    for (int y = 0; y < c.ySize(); ++y) {
      for (int x = 0; x < c.xSize(); ++x) {
        c(x, y) = (c_tmp(2*x, 2*y) + c_tmp(2*x+1, 2*y) + c_tmp(2*x, 2*y+1) + c_tmp(2*x+1, 2*y+1)) / 4.0;
        d_h(x, y) = (c_tmp(2*x, 2*y) - c_tmp(2*x+1, 2*y) + c_tmp(2*x, 2*y+1) - c_tmp(2*x+1, 2*y+1)) / 4.0;
        d_v(x, y) = (c_tmp(2*x, 2*y) + c_tmp(2*x+1, 2*y) - c_tmp(2*x, 2*y+1) - c_tmp(2*x+1, 2*y+1)) / 4.0;
        d_d(x, y) = (c_tmp(2*x, 2*y) - c_tmp(2*x+1, 2*y) - c_tmp(2*x, 2*y+1) + c_tmp(2*x+1, 2*y+1)) / 4.0;
      }
    }

    c_tmp = c;

    if (borders) {
      drawBorder(c, 2);
      drawBorder(d_h, 2);
      drawBorder(d_v, 2);
      drawBorder(d_d, 2);
    }

    output.paste(c, 0, 0);
    output.paste(d_h, c.xSize(), 0);
    output.paste(d_v, 0, c.ySize());
    output.paste(d_d, c.xSize(), c.ySize());
  }
}

// do reverse Haar decomposition
void recompHaar(CMatrix<float> & input, CMatrix<float> & output, const int & lvl) {
  assert(input.xSize() == output.xSize());
  assert(input.ySize() == output.ySize());

  int currX;
  int currY;
  output = input;
  CMatrix<float> c;
  CMatrix<float> d_h;
  CMatrix<float> d_v;
  CMatrix<float> d_d;

  for (int l = 1; l <= lvl; ++l) {
    currX = input.xSize()/pow(2, lvl-l);
    currY = input.ySize()/pow(2, lvl-l);

    output.cut(c, 0, 0, currX/2-1, currY/2-1);
    output.cut(d_h, currX/2, 0, currX-1, currY/2-1);
    output.cut(d_v, 0, currY/2, currX/2-1, currY-1);
    output.cut(d_d, currX/2, currY/2, currX-1, currY-1);

    for (int y = 0; y < currX/2; ++y) {
      for (int x = 0; x < currY/2; ++x) {
        output(2*x, 2*y) = (c(x, y) + d_h(x, y) + d_v(x, y) + d_d(x, y)) / 1.0;
        output(2*x+1, 2*y) = (c(x, y) - d_h(x, y) + d_v(x, y) - d_d(x, y)) / 1.0;
        output(2*x, 2*y+1) = (c(x, y) + d_h(x, y) - d_v(x, y) - d_d(x, y)) / 1.0;
        output(2*x+1, 2*y+1) = (c(x, y) - d_h(x, y) - d_v(x, y) + d_d(x, y)) / 1.0;
      }
    }
  }
}

// do wavelet shrinkage
void waveletShrinkage(const CMatrix<float> & input, CMatrix<float> & output, const int & threshold = 30, const bool & soft = true, const int & lvl = 3) {
  assert(input.xSize() == output.xSize());
  assert(input.ySize() == output.ySize());

  CMatrix<float> decomp(input);
  decompHaar(input, decomp, lvl);

  for (int y = decomp.ySize()/pow(2, lvl); y < decomp.ySize(); ++y) {
    for (int x = decomp.xSize()/pow(2, lvl); x < decomp.xSize(); ++x) {
      if (!soft) {
        // hard shrinkage
        if (abs(decomp(x, y)) < threshold)
          decomp(x, y) = 0;
      } else {
        // soft shrinkage
        if (abs(decomp(x, y)) > threshold)
          decomp(x, y) = (decomp(x, y)/abs(decomp(x, y)))*(abs(decomp(x, y)) - threshold);
        else
          decomp(x, y) = 0;
      }
    }
  }

  recompHaar(decomp, output, lvl);
}

double PSNR(const CMatrix<float> & orig, const CMatrix<float> & noise) {
  assert(orig.xSize() == noise.xSize());
  assert(orig.ySize() == noise.ySize());

  double mse = 0.0;

  for (int y = 0; y < orig.ySize(); ++y) {
    for (int x = 0; x < orig.xSize(); ++x) {
      mse += pow(orig(x, y) - noise(x, y), 2);
    }
  }

  mse /= orig.xSize()*orig.ySize();

  return 20 * log10(255) - 10 * log10(mse);
}

int main(int argc, const char* argv[]) {
  int ws_thresh = 30;
  bool use_soft = true;
  if (argc >= 2)
    ws_thresh = atoi(argv[1]);
  if (argc >= 3)
    use_soft = atoi(argv[2]);

  CMatrix<float> aImage;
  CMatrix<float> bImage;
  CMatrix<float> cImage;
  aImage.readFromPGM("Barbara.pgm");
  bImage.readFromPGM("BarbaraNoisy.pgm");
  cImage.readFromPGM("birds.pgm");
  CMatrix<float> aResult(aImage.xSize(), aImage.ySize(), 255);
  CMatrix<float> bResult(bImage.xSize(), bImage.ySize(), 255);
  CMatrix<float> cResult(cImage.xSize(), cImage.ySize(), 255);
  CMatrix<float> wsResult(bImage.xSize(), bImage.ySize(), 255);

  decompHaar(aImage, aResult, 3, true);
  decompHaar(bImage, bResult, 3, true);
  decompHaar(cImage, cResult, 3, true);

  waveletShrinkage(bImage, wsResult, ws_thresh, use_soft);
  wsResult.normalize(0, 255);

  std::cout << "PSNR original -- noise:" << std::endl;
  std::cout << PSNR(aImage, bImage) << std::endl;
  std::cout << "PSNR noise -- wavelet shrinkage:" << std::endl;
  std::cout << PSNR(bImage, wsResult) << std::endl;
  std::cout << "PSNR original -- wavelet shrinkage:" << std::endl;
  std::cout << PSNR(aImage, wsResult) << std::endl << std::endl;

  std::cout << "wavelet shrinkage: " << (use_soft ? "soft" : "hard") << std::endl;
  std::cout << "threshold: " << ws_thresh << std::endl;

  aResult.writeToPGM("BarbaraDecomp.pgm");
  bResult.writeToPGM("BarbaraNoisyDecomp.pgm");
  cResult.writeToPGM("birdsDecomp.pgm");
  wsResult.writeToPGM("BarbaraNoisyWS.pgm");
  return 0;
}
