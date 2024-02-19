// Simpix
// This code utilizes the root framework as well as the physics of the annealing process (in which a hot material is 
// slowly cooled in order to reach its minimum internal energy) to use the pixels of one image to recreate another.
// By making the "temperature" extremely hot to start off and moving the pixels around somewhat randomly, the color 
// difference between corresponding pixels on both images is analogous to the energy of particles of a material.

#include "TROOT.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TASImage.h"
#include "TApplication.h"
#include "TSystem.h"


#include "assert.h"

#include <math.h>
#include <iostream>
#include <stdio.h>
#include <iomanip>
#include <string>
#include <sstream>
using namespace std;

const UInt_t six16 = 16*16*16*16*16*16;
const UInt_t four16 = 16*16*16*16;
const UInt_t two16 = 16*16;


// Pixel object
struct Pixel {
  UInt_t r, g, b;
};


// Separation of red, blue, and green from hexcode
void sepRGBPixel(UInt_t p, Pixel *rgb) {
  UInt_t alpha;
  alpha = p/six16;
  rgb->r = (p - (alpha*six16))/four16;
  rgb->g = (p - (alpha*six16) - (rgb->r*four16))/two16;
  rgb->b = (p - (alpha*six16) - (rgb->r*four16) - (rgb->g*two16));
}


// Separates red, blue, and green for every pixel in an array
void sepRGBArray(UInt_t *arr, Pixel *arrRGB, int npixels) {
  int ip;
  for (ip = 0; ip < npixels; ip++) {
    sepRGBPixel(arr[ip], &arrRGB[ip]);
  }
}


// Find difference between 2 colors, ensuring output is positive
double energyC(UInt_t out, UInt_t tgt) {
  if (out > tgt) {
    return (double)(out - tgt);
  }
  else {
    return (double)(tgt - out);
  }
}


// Finds total color difference between 2 pixels; multiple methods included but only 1 is used
double energyP(Pixel out, Pixel tgt) {
  double E = 0;
  double rmean = ((double)(out.r+tgt.r))/2.0;
  double dr = energyC(out.r, tgt.r);
  double dg = energyC(out.g, tgt.g);
  double db = energyC(out.b, tgt.b);
  E = (2 + rmean/256.0)*dr*dr + 4*dg*dg + (2 + (255.0-rmean)/256.0)*db*db;  // This was found on https://en.wikipedia.org/wiki/Color_difference
  // ---Below are alternative methods of calculating color difference. Since this is calculated thousand and millions
  //    of times, the complexity of this calculation is important for the runtime of the program---
  /*if (rmean < 128) E = 2*energyC(out.r, tgt.r)*energyC(out.r, tgt.r) + 4*energyC(out.g, tgt.g)*energyC(out.g, tgt.g) + 3*energyC(out.b, tgt.b)*energyC(out.b, tgt.b);
  else E = 3*energyC(out.r, tgt.r)*energyC(out.r, tgt.r) + 4*energyC(out.g, tgt.g)*energyC(out.g, tgt.g) + 2*energyC(out.b, tgt.b)*energyC(out.b, tgt.b);*/
  return sqrt(E);
}


// Sum of color differences of all pixels
double energyT(Pixel *out, Pixel *tgt, int npixels) {
  double E = 0;
  for (int i = 0; i < npixels; i++) {
    E += energyP(out[i], tgt[i]);
  }
  return E;
}


// Swaps the reds, greens, or blues of 2 pixels
void switchColors(UInt_t *c1, UInt_t *c2) {
  UInt_t temp = *c1;
  *c1 = *c2;
  *c2 = temp;
}


// Swaps the colors of 2 pixels
void switchPixels(Pixel *p1, Pixel *p2) {
  switchColors(&p1->r, &p2->r);
  switchColors(&p1->g, &p2->g);
  switchColors(&p1->b, &p2->b);
}


// Initialize by simulating a "high temperature state" where pixels swap randomly and frequently
void InitializeHot(Pixel *out, int npixels) { 
  int i;
  for (i = npixels - 1; i > 0; i--) {
    int j = ((int)(npixels*drand48())) % (i+1);
    switchPixels(&out[i], &out[j]);
  }
}


// Use Boltzmann distribution as probability to swap 2 pixels according to the "temperature"; Used 50% of the time
void updatePixels(Pixel *out, Pixel *tgt, int npixels, double beta, int n1, int n2, double &E) {
  double dE = energyP(out[n2], tgt[n1]) + energyP(out[n1], tgt[n2]) - energyP(out[n1], tgt[n1]) - energyP(out[n2], tgt[n2]);
  if (dE < 0 || drand48() < exp(-dE*beta)) {
    switchPixels(&out[n1], &out[n2]);
    E += dE;
  }
}


// Use Boltzmann distribution as probability to reverse subsequence of pixels according to the "temperature"; Used 50% of the time
void updatePixels2(Pixel *out, Pixel *tgt, int npixels, double beta, int n1, int n2, double &E) {
  int diff = abs(n2-n1) + 1;
  double dE = 0;
  for (int i = 0; i < diff/2; i++) {
    dE += energyP(out[n2-i], tgt[n1+i]) + energyP(out[n1+i], tgt[n2-i]) - energyP(out[n1+i], tgt[n1+i]) - energyP(out[n2-i], tgt[n2-i]);
  }
  if (dE < 0 || drand48() < exp(-dE*beta)) {
    for (int i = 0; i < diff/2; i++) {
      switchPixels(&out[n1+i], &out[n2-i]);
    }
    E += dE;
  }
}


// Process of updating the "state" or image based on "energy" or total color difference and "temperature"
void sweep(Pixel *out, Pixel *tgt, int npixels, double beta, double &E) {
  int np, i, j;
  for (np = 0; np < npixels; np++) {
    if (drand48() > 0.5) {
      i = (int)(drand48()*npixels);
      j = (int)(drand48()*npixels);
      while (i == j) {
        j = (int)(drand48()*npixels);
      }
      updatePixels(out, tgt, npixels, beta, i, j, E);
    } else {
      i = (int)(drand48()*npixels);
      if (i < npixels-300) {
        j = i + (int)(drand48()*200);
        while (j == i) { 
          j = i + (int)(drand48()*200);
        }
      } else if (i < npixels-50) {
        j = i + (int)(drand48()*40);
        while (j == i) {
          j = i + (int)(drand48()*40);
        }
      } else {
        while (i > npixels-20) {
          i -= drand48()*20;
        }
        j = i + (int)(drand48()*20);
        while (j == i) {
          j = i + (int)(drand48()*20);
        }
      }
      updatePixels2(out, tgt, npixels, beta, i, j, E);
    }
  }
}


int main(int argc, char **argv){

  if (argc<3) {
    cout << "Usage: simpix image1 image2 <output=out.png>" << endl;
    return 0; 
  }
  TString fsrc=argv[1];   // This is the image which will have its pixels manipulated
  fsrc = fsrc;
  TString ftgt=argv[2];   // This is the image which will be copied
  ftgt = ftgt;
  TString fout;
  argc>3 ? fout = argv[3] : fout="out.png";
  fout = "Output/" + fout;
  cout << "Reading images: source= " << fsrc << " target= " << ftgt << endl;
  cout << "Output= " << fout << endl;

  TString energyFileName = "energy.dat";
  string collageFileEnding = "collage.png";

  // Create image objects
  TASImage *src = new TASImage(fsrc.Data());
  TASImage *tgt = new TASImage(ftgt.Data());
  TASImage *out = new TASImage(*src); // start with copy of source

  // Test image geometry, exit if they are not the same dimensions
  assert ( src->GetWidth() == tgt->GetWidth() && src->GetHeight() == tgt->GetHeight() );
  cout << "Pixel Geometry: " << src->GetWidth() << " x " << src->GetHeight() << endl;
  cout << "Total Pixels: " << (long)src->GetWidth()*(long)src->GetHeight() << endl;
  Long_t numPix=(long)src->GetWidth()*(long)src->GetHeight();

  // *** The work happens here
  // access the pixels for the output image 
  // each pixel is a 32-bit word, 1 byte each for (alpha,red,green,blue)
  // don't touch alpha (bits 31:28)
  UInt_t *outPix = out->GetArgbArray();
  UInt_t *tgtPix = tgt->GetArgbArray();
  // Previous issues with running out of memory for larger images, which dynamic memory allocation seems to have fixed
  Pixel* outRGB = new Pixel[numPix];
  Pixel* tgtRGB = new Pixel[numPix];
  sepRGBArray(outPix, outRGB, numPix);
  sepRGBArray(tgtPix, tgtRGB, numPix);

  InitializeHot(outRGB, numPix);    // Start the "state" at a high "temperature" so that random fluctuations are extremely likely

  // *************************

  // Begin "annealing" process of the image
  // it = iterative variable in loop of decreasing "temperature"; nt = total # of "temperature" steps for "annealing" process
  // itherm = iterative variable for updating the image; ntherm = total # of times the image is updated
  int it, nt, itherm, ntherm;
  double T, beta, Tmax, E;
  nt = 100;
  ntherm = 10;
  Tmax = 1000;

  FILE *output;
  output = fopen(energyFileName, "w");

  E = energyT(outRGB, tgtRGB, numPix);
  for (it = nt; it > 0; it--) {   // "Temperature" decrement loop
    T = (Tmax*((double)it))/((double)nt);
    beta = 1/T;
    for (itherm = 0; itherm < ntherm; itherm++) { sweep(outRGB, tgtRGB, numPix, beta, E); }// Allow the "state" to "thermalize" at this "temperature"
    fprintf(output, "%lf  %lf\n", T, E);    // Print "temperature" and "energy" values to the output file
    double progress = (double)(nt-it)/(double)nt*100.;  // These last few lines of the loop are for printing progress
    cout << "\rProgress: " << std::fixed << setprecision(2) << progress << "%";
    cout.flush();
  }
  cout << endl;
  for (int i = 0; i < numPix; i++) {    // Convert final pixels back to hexcode
    UInt_t hexcode = 0;
    hexcode += 255*six16;
    hexcode += outRGB[i].r*four16;
    hexcode += outRGB[i].g*two16;
    hexcode += outRGB[i].b;
    outPix[i] = hexcode;
  }



  // print the results
  TCanvas *canvas1 = new TCanvas("canvas1", "images", 640, 480);
  canvas1->Divide(2,2);

  canvas1->cd(1);
  canvas1->Draw();
  src->Draw("X");
  canvas1->cd(2);
  tgt->Draw("X");
  canvas1->cd(3);
  out->Draw("X");
  stringstream ssw, ssh;
  ssw << src->GetWidth();
  ssh << src->GetHeight();
  string collageFileName = "Output/"+ssw.str()+"x"+ssh.str()+collageFileEnding;
  canvas1->Print(collageFileName.c_str());
  
  // save the new image
  out->WriteImage(fout.Data());

  // comment out the lines for running in batch mode
  // cout << "Press ^c to exit" << endl;
  // theApp.SetIdleTimer(30,".q");  // set up a failsafe timer to end the program  
  // theApp.Run();

  delete outRGB;
  delete tgtRGB;
}
