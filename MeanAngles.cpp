// MeanAngles.cpp : Diese Datei enthält die Funktion "main". Hier beginnt und endet die Ausführung des Programms.
//

//Including root functionalities:
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TMultiGraph.h>
#include <TLine.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TStyle.h>
#include <TString.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TMath.h>
#include <TLegend.h>
#include <TCut.h>
#include <THStack.h>
#include <TGaxis.h>
#include <TF1.h>
#include <TError.h> // root verbosity level
#include <TApplication.h>
#include <TNtuple.h>
#include <TImage.h>
#include <TAttImage.h>
#include <TPaveLabel.h>
#include <TLatex.h>

//Including standard C/C++-Libraries:
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/types.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <sstream>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <ctime>

using namespace std;

int main(int argc, char *argv[])
{
    const int maxFileNumber = 100;  //Arbitrary number; adapt if more than 100 files need to be analyzed
    const int nCh = 8;
    const int nBeamPositions = 12;
    const int nBeamPositions30 = 24;
    cout << endl;
    string outPath = "../runs/MeanAngles/";
    string inPath = "../runs/PhotonCountDistribution/";
    TTree* tree[maxFileNumber];
    int angleList[] = { 135,180,225,270,90,45,0,315 }; //Represents channels 0-7
    float womX = 310;   //Position of WOM D
    float womY = -510;

    //For 0 degrees:
    float positionX[nBeamPositions] = { 0,160,320,-320,-160,-160,-320,320,160,160,204,310 };        //Corresponding to Positions 0 - 11;
    float positionY[nBeamPositions] = { 0,160,320,320,160,-160,-320,-320,-160,-510,-404,-360 };
    float womDistance[nBeamPositions];  //Distances of positions 0 - 11 to the center of WOM D (in mm)
    float womAngle[nBeamPositions];     //Angles of positions relative to the center of WOM D (in degrees)

    //For 30 degrees:
    float positionX30[nBeamPositions30];
    float drawnPositionX30[nBeamPositions30] = { 0,140,282,-282,-140,-140,-282,282,140,160,204,310,-282,282,0,0,0,0,-140,140,282,240,0,-290 };   //Corresponding to positions 0 - 23;
    //Calculating the "average beam position" (adjusting for the 30° tilt, see bachelor thesis):
    for (int i = 0; i < nBeamPositions30; i++)
    {
        positionX30[i] = drawnPositionX30[i] + 96.7; 
    }
    float positionY30[nBeamPositions30] = { 0,160,320,320,160,-160,-320,-320,-160,-510,-404,-360,0,0,407,-460,-160,160,0,0,-510,-510,-510,-510 };
    float womDistance30[nBeamPositions30];
    float womAngle30[nBeamPositions30];

    //Computing angles and distances from beam positions to WOM D center:
    for (int i = 0; i < nBeamPositions; i++)
    {
        womDistance[i] = TMath::Sqrt(TMath::Power(positionX[i] - womX, 2) + TMath::Power(positionY[i] - womY, 2));  
        if (positionX[i] >= womX)
        {
            womAngle[i] = (TMath::ATan((positionY[i] - womY) / (positionX[i] - womX))) * 180.0 / TMath::Pi();
        }
        else
        {
            womAngle[i] = (TMath::ATan((positionY[i] - womY) / (positionX[i] - womX)) + TMath::Pi()) * 180.0 / TMath::Pi();
        }
    }
    for(int i = 0; i < nBeamPositions30; i++)
    {
        womDistance30[i] = TMath::Sqrt(TMath::Power(positionX30[i] - womX, 2) + TMath::Power(positionY30[i] - womY, 2));
        if (positionX30[i] >= womX)
        {
            womAngle30[i] = (TMath::ATan((positionY30[i] - womY) / (positionX30[i] - womX))) * 180.0 / TMath::Pi();
        }
        else
        {
            womAngle30[i] = (TMath::ATan((positionY30[i] - womY) / (positionX30[i] - womX)) + TMath::Pi()) * 180.0 / TMath::Pi();
        }
    }

    //Reading total photon counts for each run from .txt-file created by PhotonCountDistribution:
    string photonCountString[maxFileNumber];
    string photonCountBuffer;
    string photonCountPath = Form("%s%s.txt", inPath.c_str(), "photonCounts");
    ifstream photonCountFile;
    photonCountFile.open(photonCountPath, ifstream::in);
    if (photonCountFile.is_open())
    {
        int argStr = 0;
        while (getline(photonCountFile, photonCountBuffer))
        {
            photonCountString[argStr] = photonCountBuffer;
            argStr += 1;
        }
        photonCountFile.close();
    }

    //Parameters, canvases and graphs comparing peak positions to beam positions:

    //For Peak Location vs Beam Angle graph:
    float peakValues[maxFileNumber][nCh + 1];
    float peakValueErrors[maxFileNumber][nCh + 1];
    float beamAngles[maxFileNumber];

    //For Standard Deviation vs Beam Angle graph:
    float stdDevs[maxFileNumber][nCh + 1];
    float stdDevErrors[maxFileNumber][nCh + 1];

    //For Peak-Valley-Ratio vs light yield Graph (formerly vs beam distances, beam distances were kept because they might be needed):
    float nPhotons[maxFileNumber];
    float beamDistances[maxFileNumber];
    float peakValleyRatio[maxFileNumber];
    float peakValleyRatioErrors[maxFileNumber];

    TLatex* latex[maxFileNumber];

    TCanvas peakCanvas("peakCanvas", "Peak Position Comparison", 1557, 2000);
    TCanvas peakCanvas8("peakCanvas8", "Peak Positions All Channels", 1920, 1080);
    peakCanvas.cd();
    TPaveLabel peakTitle(0.1, 0.96, 0.9, 0.99, "Comparison of #phi_{peak} values across runs");
    peakTitle.SetTextSize(.7);
    peakTitle.SetLineColor(0);
    peakTitle.SetBorderSize(0);
    peakTitle.SetFillColor(0);
    peakTitle.Draw();
    TPad peakPad("Graphs", "Graphs", 0, 0, 1, 0.96);
    peakPad.Draw();
    peakPad.cd();
    peakPad.Divide(2, 4);
    //peakPad.SetLeftMargin(.2);
    //peakPad.SetBottomMargin(.15);

    TCanvas peakValleyCanvas0("peakValleyCanvas0", "Peak-Valley-Ratio Comparison 0 Degrees", 1920, 1080);
    TCanvas peakValleyCanvas30("peakValleyCanvas0", "Peak-Valley-Ratio Comparison 30 Degrees", 1920, 1080);
    TCanvas sigmaCanvas0("sigmaCanvas0", "Sigma Comparison for 0 Degrees", 1920, 1080);
    TCanvas sigmaCanvas30("sigmaCanvas0", "Sigma Comparison for 30 Degrees", 1920, 1080);
    TCanvas correlationCanvas("correlation Canvas", "Correlation between sigma and r values", 1920, 1080);
    TCanvas canvas91011("canvas91011", "Comparison between Positions 9, 10 and 11 at 1.4 GeV", 1920, 1080);

    TCanvas stdDevCanvas("stdDevCanvas", "Standard Deviation Comparison", 1557, 2000);
    TCanvas stdDevCanvas8("stdDevCanvas8", "Standard Deviations All Channels", 1920, 1080);
    stdDevCanvas.cd();
    TPaveLabel stdDevTitle(0.1, 0.96, 0.9, 0.99, "Comparison of standard deviations across runs");
    stdDevTitle.SetTextSize(.7);
    stdDevTitle.SetLineColor(0);
    stdDevTitle.SetBorderSize(0);
    stdDevTitle.SetFillColor(0);
    stdDevTitle.Draw();
    TPad stdDevPad("Graphs", "Graphs", 0, 0, 1, 0.96);
    stdDevPad.Draw();
    stdDevPad.cd();
    stdDevPad.Divide(2, 4);
    //stdDevPad.SetLeftMargin(.2);
    //stdDevPad.SetBottomMargin(.15);

    //For peak locations, two graphs are necessary for 0 and 30 degrees:
    TGraphErrors *peakGraph0[nCh + 1];
    TGraphErrors *peakGraph30[nCh + 1];
    TMultiGraph *peakMultiGraph[nCh + 1];
    TLegend *peakLegend[nCh + 1];
    float peakLegendXMin[nCh + 1] = { 0.535,0.535,0.535,0.535,0.535,0.535,0.535,0.535,0.5 };
    float peakLegendXMax[nCh + 1] = { 0.935,0.935,0.935,0.935,0.935,0.935,0.935,0.935,0.9 };
    float peakLegendYMin[nCh + 1] = { 0.65,0.65,0.052,0.052,0.65,0.65,0.052,0.052,0.38 };
    float peakLegendYMax[nCh + 1] = { 0.9,0.9,0.302,0.302,0.9,0.9,0.302,0.302,0.63 };

    //One Graph for Correlation between sigma and r values:
    TGraphErrors *correlationGraph0 = new TGraphErrors;
    TGraphErrors *correlationGraph30 = new TGraphErrors;
    TMultiGraph* correlationMultiGraph = new TMultiGraph("correlationMultiGraph", "Correlation between #it{r} and #sigma values for individual runs;#sigma;#it{r}");
    TLegend *correlationLegend = new TLegend(0.5, 0.6, 0.9, 0.9);

    //Six different Graphs are created and merged for the peak Valley ratios (for three different energies and two angles):
    TGraphAsymmErrors *peakValleyGraph0_14 = new TGraphAsymmErrors();
    TGraphAsymmErrors *peakValleyGraph0_26 = new TGraphAsymmErrors();
    TGraphAsymmErrors *peakValleyGraph0_52 = new TGraphAsymmErrors();
    TGraphAsymmErrors *peakValleyGraph30_14 = new TGraphAsymmErrors();
    TGraphAsymmErrors *peakValleyGraph30_26 = new TGraphAsymmErrors();
    TGraphAsymmErrors *peakValleyGraph30_52 = new TGraphAsymmErrors();
    TMultiGraph *peakValleyMultiGraph0 = new TMultiGraph();
    TMultiGraph *peakValleyMultiGraph30 = new TMultiGraph();
    TLegend *peakValleyLegend0 = new TLegend(0.18,0.6,0.58,0.9);
    TLegend *peakValleyLegend30 = new TLegend(0.18, 0.6, 0.58, 0.9);

    TGraphErrors* sigmaGraph0_14 = new TGraphErrors();
    TGraphErrors* sigmaGraph0_26 = new TGraphErrors();
    TGraphErrors* sigmaGraph0_52 = new TGraphErrors();
    TGraphErrors* sigmaGraph30_14 = new TGraphErrors();
    TGraphErrors* sigmaGraph30_26 = new TGraphErrors();
    TGraphErrors* sigmaGraph30_52 = new TGraphErrors();
    TMultiGraph* sigmaMultiGraph0 = new TMultiGraph();
    TMultiGraph* sigmaMultiGraph30 = new TMultiGraph();
    TLegend* sigmaLegend0 = new TLegend(0.5, 0.6, 0.9, 0.9);
    TLegend* sigmaLegend30 = new TLegend(0.5, 0.6, 0.9, 0.9);
    
    //Two separate graphs for stdDev are combined into a multigraph, so that points for 0 and 30 degrees can be differently coloured:
    TGraphErrors *stdDevGraph0[nCh + 1];
    TGraphErrors *stdDevGraph30[nCh + 1];
    TMultiGraph *stdDevMultiGraph[nCh + 1];
    TLegend *stdDevLegend[nCh + 1];

    //For comparison of positions 9, 10 and 11:
    TLegend* legend91011 = new TLegend(0.6, 0.7, 0.9, 0.9);
    TLine* line9 = new TLine(0, 0, 1, 1);
    TLine* line10 = new TLine(0, 0, 1, 1);
    TLine* line11 = new TLine(0, 0, 1, 1);
    line9->SetLineColor(433);
    line10->SetLineColor(617);
    line11->SetLineColor(797);
    legend91011->AddEntry(line9, "Position 9", "l");
    legend91011->AddEntry(line10, "Position 10", "l");
    legend91011->AddEntry(line11, "Position 11", "l");

    //Style Settings:
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetGridStyle(3);
    gStyle->SetGridWidth(1);
    gStyle->SetGridColor(16);
    gStyle->SetLineScalePS(1);

    //Creating variables to read data from the TTree stored in the .root file into our own TTree:
    int runNumber[maxFileNumber];
    int runPosition[maxFileNumber];
    float runEnergy[maxFileNumber];
    int runAngle[maxFileNumber];
    float photonCount[maxFileNumber];

    //The .root files are passed to this program as a list of arguments by the corresponding .sh script. This is the beginning of the loop over the .root files:
    for (int arg = 1; arg < argc; arg++)
    {
        nPhotons[arg] = stof(photonCountString[arg - 1].c_str());
        tree[arg] = new TTree;
        string filePath = argv[arg];
        printf("Analyzing file %s\n", filePath.c_str());

        //Open .root file:
        TFile file(filePath.c_str());
        if (file.IsZombie())
        {
            cout << "Problem with file " << filePath << "; check if file path is correct!" << endl;
            exit(-1);
        }

        //Reading Data from the tree named "T" in the root file into our own TTree:
        file.GetObject("T", tree[arg]);
        tree[arg]->SetBranchAddress("runNumber", &runNumber[arg]);
        tree[arg]->SetBranchAddress("runPosition", &runPosition[arg]);
        tree[arg]->SetBranchAddress("runEnergy", &runEnergy[arg]);
        tree[arg]->SetBranchAddress("runAngle", &runAngle[arg]);
        tree[arg]->GetEntry(1);

        //Storing information for later use when comparing beam angles/distances and peak positions:
        if (runAngle[arg] == 0)
        {
            beamAngles[arg] = womAngle[runPosition[arg]];
            beamDistances[arg] = womDistance[runPosition[arg]];

        }
        else
        {
            beamAngles[arg] = womAngle30[runPosition[arg]];
            beamDistances[arg] = womDistance30[runPosition[arg]];
        }

        //Reading the photonCut threshold calculated by PhotonCountDistribution.cpp:
        /*
        string photonCut[nCh];
        string buffer;
        string cutoffPath = Form("%s%s%d.txt", inPath.c_str(), "CutoffValues_run", runNumber[arg]);
        ifstream cutoffFile;
        cutoffFile.open(cutoffPath, ifstream::in);
        if (cutoffFile.is_open())
        {
            int nString = 0;
            while (getline(cutoffFile, buffer))
            {
                photonCut[nString] = buffer;
                nString += 1;
            }
            cutoffFile.close();
        }
        */
        //string selection = Form("Integral[0] > %s && Integral[1] > %s && Integral[2] > %s && Integral[3] > %s && Integral[4] > %s && Integral[5] > %s && Integral[6] > %s && Integral[7] > %s && chargeChannelSumWOM[3] > 5",
            //photonCut[0].c_str(), photonCut[1].c_str(), photonCut[2].c_str(), photonCut[3].c_str(), photonCut[4].c_str(), photonCut[5].c_str(), photonCut[6].c_str(), photonCut[7].c_str());

        string selection = "chargeChannelSumWOM[3] > 10";
        //Canvas for 9 histograms:
        TCanvas canvas;
        TCanvas canvasAllCh;
        
        TH1F* histVec[9];
        int histCtr = 0;
        TLegend* legend[9];
        TF1* fitFunction[9];
        int lowerHistEdge[9] = { -225,-180,-135,-90,-270,-315,0,-45,-180 };
        int upperHistEdge[9] = { 134,179,224,269,89,44,359,314,179 };
        int lowerFitEdge, upperFitEdge;

        //Loop for drawing Histograms -- draws the Histogram corresponding to weightedMeanAngles[i], see read.C line >~1000
        canvasAllCh.SetCanvasSize(1920, 1080);
        canvasAllCh.SetTitle(Form("Run %d - All Channels", runNumber[arg]));
        canvasAllCh.SetGrid();
        canvasAllCh.cd();
        gPad->SetLeftMargin(.18);
        gPad->SetBottomMargin(.16);
        canvas91011.cd();
        gPad->SetLeftMargin(.18);
        gPad->SetBottomMargin(.18);
        canvas.cd();
        canvas.SetCanvasSize(1577, 2000);
        canvas.SetTitle(Form("Run %d", runNumber[arg]));
        canvas.cd(0);
        TPaveLabel title(0.1, 0.96, 0.9, 0.99, Form("Distributions of #bar{#phi_{ew}}, Run %d, Pos. %d, %1.1f GeV", runNumber[arg], runPosition[arg], runEnergy[arg] / 10.0));
        TPaveLabel xTitle(0, 0.01, 1, 0.03, "Event-Wise Mean Angle #bar{#phi_{ew}} [deg.]");
        TPaveLabel yTitle(0.01, 0, 0.03, 1, "Number of Entries");
        title.SetTextSize(.7);
        xTitle.SetTextSize(.7);
        yTitle.SetTextAngle(90);
        yTitle.SetTextSize(.017);
        title.SetLineColor(0);
        xTitle.SetLineColor(0);
        yTitle.SetLineColor(0);
        title.SetBorderSize(0);
        xTitle.SetBorderSize(0);
        yTitle.SetBorderSize(0);
        title.SetFillColor(0);
        xTitle.SetFillColor(0);
        yTitle.SetFillColor(0);
        title.Draw();
        xTitle.Draw();
        yTitle.Draw();
        TPad pad("Graphs", "Graphs", 0.03, 0.03, 1, 0.96);
        pad.Draw();
        pad.cd();
        pad.Divide(2, 4);
        //pad.SetLeftMargin(0.065);
        //pad.SetBottomMargin(0.15);

        for (int i = 0; i < 9; i++)
        {
            TString histName, histTitle, histDraw;
            //Drawing the histograms:
            if (i == 8)
            {
                histTitle.Form("Distribution of #bar{#phi_{ew}} for all Channels, Run %d, Pos. %d, %1.1f GeV;#bar{#phi_{ew}} [deg.];Number of Entries", runNumber[arg], runPosition[arg], runEnergy[arg] / 10.0);
                canvasAllCh.cd();
                legend[i] = new TLegend(0.55, 0.6, 0.9, 0.9);
                //legend[i]->AddEntry(histVec[i], "Weighted mean Angles");
                gPad->SetBottomMargin(.17);
            }
            else
            {
                histTitle.Form("Ignoring Channel %d (%d deg.)", i, angleList[i]); //;#bar{#phi_{ew}} [deg.];Number of Entries
                pad.cd(i + 1);
                legend[i] = new TLegend(0.5585, 0.6, 0.935, 0.9);
                gPad->SetLeftMargin(.065);
                gPad->SetBottomMargin(.052);
                gPad->SetRightMargin(.065);
                gPad->SetGrid();
                TString histName, histTitle, histDraw;
            }
            histName.Form("Hist%d", i);
            histDraw.Form("weightedMeanAngle[%d]>>Hist%d", i, i);
            histVec[i] = new TH1F(histName, histTitle, 360, lowerHistEdge[i], upperHistEdge[i]);
            tree[arg]->Draw(histDraw, selection.c_str(), "HISTE");
            legend[i]->SetFillColorAlpha(kWhite, .85);
            //legend[i]->AddEntry(histVec[i], "Weighted mean Angles", "lpf");
            //Fill the histogram with the entries from the weightedMeanAngle branch, but only with entries where ALL channels contain photon counts larger than the threshold calculated in
            //PhotonCountDistribution.cpp.

            //Calculating the fit functions for all channels included (combination of 3 gaussians) and the ignored channels (gaussians):
            lowerFitEdge = lowerHistEdge[i];
            upperFitEdge = upperHistEdge[i];
            if (i == 8)
            {
                canvasAllCh.cd();
                histVec[i]->SetLineWidth(3.5);
                //Getting estimators for starting values:
                float peakStartingValue = histVec[i]->GetMean();
                float stdDevStartingValue = histVec[i]->GetStdDev();
                float offsetStartingValue = histVec[i]->GetMinimum();
                float scalingStartingValue = histVec[i]->GetMaximum();

                //Creating the fit function, setting the starting parameters to these estimators and fitting it to the histogram:
                fitFunction[i] = new TF1("TripleGaussFit", "[0]*(exp(-0.5*((x-[1]+360.0)/[2])^2)+exp(-0.5*((x-[1]-360.0)/[2])^2)+exp(-0.5*((x-[1])/[2])^2))+[3]", -540, 540);
                fitFunction[i]->SetParameters(scalingStartingValue, peakStartingValue, stdDevStartingValue, offsetStartingValue);
                fitFunction[i]->SetParLimits(0, 0, 10000);   //Making sure the gauss function is positive
                fitFunction[i]->SetParLimits(1, -180, 180); //Making sure the central fit peak is in the histogram range
                fitFunction[i]->SetParLimits(2, 0, 10000);   //Making sure the standard deviation is positive (could also just use absolute value); upper limits are arbitrary but should cover everything
                fitFunction[i]->SetLineWidth(5);
                histVec[i]->Fit("TripleGaussFit", "RQM");

                histVec[i]->GetXaxis()->SetTitleSize(.07);
                histVec[i]->GetXaxis()->SetLabelSize(.06);
                histVec[i]->GetYaxis()->SetTitleSize(.07);
                histVec[i]->GetYaxis()->SetLabelSize(.06);
                
                //Drawing the histogram with the corresponding part of the fit function:
                histVec[i]->DrawCopy();
                fitFunction[i]->Draw("same");

                if (runNumber[arg] == 31)
                {
                    canvas91011.cd(0);
                    histVec[i]->SetLineColor(433);
                    histVec[i]->SetTitle("Comparison between Positions 9, 10 and 11 for all Channels at 1.4 GeV;#bar{#phi_{ew}} [deg.];Number of Entries");
                    histVec[i]->DrawCopy();

                }
                else if (runNumber[arg] == 35)
                {
                    canvas91011.cd(0);
                    histVec[i]->SetLineColor(617);
                    histVec[i]->DrawCopy("same");
                    
                }
                else if (runNumber[arg] == 36)
                {
                    canvas91011.cd(0);
                    histVec[i]->SetLineColor(797);
                    histVec[i]->DrawCopy("same");
                }

                canvasAllCh.cd();
                gPad->SetLeftMargin(.18);
                gPad->SetBottomMargin(.16);
                //Extracting data from the fit function:
                
                //Chi^2 / ndf:
                float rChi2 = fitFunction[i]->GetChisquare() / fitFunction[i]->GetNDF();

                //Peak Positions and errors:
                peakValues[arg][i] = fitFunction[i]->GetParameter(1);
                peakValueErrors[arg][i] = fitFunction[i]->GetParError(1);

                //Ratio of the y values of the central fit maximum and the minimum between the central and the left maximum:
                float fitMax = fitFunction[i]->GetMaximum(fitFunction[i]->GetParameter(1) - 1, fitFunction[i]->GetParameter(1) + 1); //Searches for maximum in the region around the central peak
                float fitMin = fitFunction[i]->Eval(fitFunction[i]->GetParameter(1) - 180.0); //Two maxima are 360 deg apart and the minimum is in the middle
                peakValleyRatio[arg] = fitMax / fitMin;

                //Gaussian error propagation (some of these values are retrieved again later, this is for better overview and to keep the error propagation "separate":
                float scaling = fitFunction[i]->GetParameter(0);
                float scalingError = fitFunction[i]->GetParError(0);
                float sigma = fitFunction[i]->GetParameter(2);
                float sigmaError = fitFunction[i]->GetParError(2);
                float offset = fitFunction[i]->GetParameter(3);
                float offsetError = fitFunction[i]->GetParError(3);
                float sigma2 = TMath::Power(sigma, 2);

                //Derivatives of f_max / f_min with respect to parameters (calculated by Mathematica, which is great for symbolic calculations but nothing else:
                float scalingDerivative = (1 + 2 * TMath::Exp(-32400.0 / sigma2)) / ((TMath::Exp(-72900.0 / sigma2) + 2 * TMath::Exp(-8100.0 / sigma2)) * scaling + offset)
                    - ((TMath::Exp(-72900.0 / sigma2) + 2 * TMath::Exp(-8100.0 / sigma2)) * ((1 + 2 * TMath::Exp(-32400.0 / sigma2)) * scaling + offset)) / TMath::Power(((TMath::Exp(-72900.0 / sigma2) + 2 * TMath::Exp(-8100.0 / sigma2)) * scaling + offset), 2);
                float sigmaDerivative = 129600.0 * TMath::Exp(-32400.0 / sigma2) * scaling / (TMath::Power(sigma, 3) * ((TMath::Exp(-72900.0 / sigma2) + 2 * TMath::Exp(-8100.0 / sigma2)) * scaling + offset))
                    - scaling * ((145800 * TMath::Exp(-72900.0 / sigma2) + 32400.0 * TMath::Exp(-8100.0 / sigma2)) / TMath::Power(sigma, 3)) * ((1 + 2 * TMath::Exp(-32400.0 / sigma2)) * scaling + offset)
                    / TMath::Power((TMath::Exp(-72900.0 / sigma2) + 2 * TMath::Exp(-8100.0 / sigma2)) * scaling + offset, 2);
                float offsetDerivative = 1 / ((TMath::Exp(-72900.0 / sigma2) + 2 * TMath::Exp(-8100.0 / sigma2)) * scaling + offset)
                    - ((1 + 2 * TMath::Exp(-32400.0 / sigma2)) * scaling + offset) / TMath::Power((TMath::Exp(-72900.0 / sigma2) + 2 * TMath::Exp(-8100 / sigma2)) * scaling + offset, 2);

                peakValleyRatioErrors[arg] = TMath::Sqrt(TMath::Power(scalingDerivative * scalingError, 2) + TMath::Power(sigmaDerivative * sigmaError, 2) + TMath::Power(offsetDerivative * offsetError, 2));

                //Standard Deviations of central fit peaks:
                stdDevs[arg][i] = fitFunction[i]->GetParameter(2);
                stdDevErrors[arg][i] = fitFunction[i]->GetParError(2);

                //Legend:
                //legend[i]->AddEntry(fitFunction[i], "Fit Function:", "l");
                legend[i]->AddEntry((TObject*)0, Form("#chi^{2} / ndf: \t%1.2f", rChi2), "");
                legend[i]->AddEntry((TObject*)0, Form("Centr. Fit Max. at: \t%1.2f #pm%1.2f", peakValues[arg][i], peakValueErrors[arg][i]), "");
                legend[i]->AddEntry((TObject*)0, Form("Centr. Std. Dev.: \t%1.2f #pm%1.2f", stdDevs[arg][i], stdDevErrors[arg][i]), "");
                legend[i]->AddEntry((TObject*)0, Form("f_{max} / f_{min}: \t%1.2f #pm%1.2f", peakValleyRatio[arg], peakValleyRatioErrors[arg]), "");
                legend[i]->AddEntry((TObject*)0, Form("Entries: %1.f", histVec[i]->GetEntries()), "");
                legend[i]->Draw("same");
                gStyle->SetTitleSize(0.08, "t");
                
            }
            else
            {
                pad.cd(i + 1);
                //Simple Gaussian Fit:
                string fitName = Form("GaussFit[%d]", i);
                fitFunction[i] = new TF1(fitName.c_str(), "gaus", lowerFitEdge, upperFitEdge);
                fitFunction[i]->SetParameters(histVec[i]->GetMaximum(), histVec[i]->GetMean(), histVec[i]->GetRMS());
                histVec[i]->Fit(fitName.c_str(), "RQM");
                histVec[i]->Draw();
                fitFunction[i]->Draw("same");
                peakValues[arg][i] = fitFunction[i]->GetParameter(1);
                peakValueErrors[arg][i] = fitFunction[i]->GetParError(1);
                stdDevs[arg][i] = fitFunction[i]->GetParameter(2);
                stdDevErrors[arg][i] = fitFunction[i]->GetParError(2);
                float rChi2 = fitFunction[i]->GetChisquare() / fitFunction[i]->GetNDF();
                //legend[i]->AddEntry(fitFunction[i], "Gaussian Fit Function", "l");
                legend[i]->AddEntry((TObject*)0, Form("Fit Mean: \t%1.2f #pm%1.2f", peakValues[arg][i], peakValueErrors[arg][i]), "");
                legend[i]->AddEntry((TObject*)0, Form("Std. Dev.: \t%1.2f #pm%1.2f", stdDevs[arg][i], stdDevErrors[arg][i]), "");
                legend[i]->AddEntry((TObject*)0, Form("#chi^{2} / ndf: \t%1.2f", rChi2), "");
                legend[i]->Draw("same");
                gStyle->SetTitleSize(0.08, "t");
                //histVec[i]->GetXaxis()->SetTitleSize(.07);
                histVec[i]->GetXaxis()->SetLabelSize(.06);
                //histVec[i]->GetYaxis()->SetTitleSize(.07);
                histVec[i]->GetYaxis()->SetLabelSize(.06);
            }
        }   //END OF LOOP OVER CHANNELS
        canvas.SaveAs(Form("%s%s%d%s%d%s%2.0f%s%d%s.pdf", outPath.c_str(), "run", runNumber[arg], "_pos", runPosition[arg], "_e", runEnergy[arg], "_", runAngle[arg], "deg"));
        canvasAllCh.SaveAs(Form("%s%s%d%s%d%s%2.0f%s%d%s_allChannels.pdf", outPath.c_str(), "run", runNumber[arg], "_pos", runPosition[arg], "_e", runEnergy[arg], "_", runAngle[arg], "deg"));
    }   //END OF LOOP OVER ARGUMENTS
    
    //Filling graphs for peak positions, peak valley ratio and standard deviations:
    int pointCounter0 = 0;
    int pointCounter0_14 = 0;
    int pointCounter0_26 = 0;
    int pointCounter0_52 = 0;
    int pointCounter30 = 0;
    int pointCounter30_14 = 0;
    int pointCounter30_26 = 0;
    int pointCounter30_52 = 0;

    for (int ch = 0; ch < nCh + 1; ch++)
    {
        //Reset the peak pointCounters for every channel:
        pointCounter0 = 0;
        pointCounter30 = 0;

        peakGraph0[ch] = new TGraphErrors();
        stdDevGraph0[ch] = new TGraphErrors();
        peakGraph30[ch] = new TGraphErrors();
        stdDevGraph30[ch] = new TGraphErrors();
        for (int arg = 1; arg < argc; arg++)
        {
            if (runAngle[arg] == 0)
            {
                if (beamAngles[arg])
                {
                    peakGraph0[ch]->SetPoint(pointCounter0, beamAngles[arg], peakValues[arg][ch]);
                    peakGraph0[ch]->SetPointError(pointCounter0, 0, peakValueErrors[arg][ch]);
                    stdDevGraph0[ch]->SetPoint(pointCounter0, beamAngles[arg], stdDevs[arg][ch]);
                    stdDevGraph0[ch]->SetPointError(pointCounter0, 0, stdDevErrors[arg][ch]);
                    pointCounter0 += 1;
                }
                if (peakValleyRatio[arg])
                {
                    if (runEnergy[arg] == 14)
                    {
                        peakValleyGraph0_14->SetPoint(pointCounter0_14, nPhotons[arg], peakValleyRatio[arg]);
                        //Manually ensuring that errors don't reach below 1, as the peakValleyRatio cannot be less than 1 per definition:
                        peakValleyGraph0_14->SetPointError(pointCounter0_14, 0, 0, TMath::Min(peakValleyRatio[arg] - 1, peakValleyRatioErrors[arg]), peakValleyRatioErrors[arg]);
                        if (nPhotons[arg] >= 1000)
                        {
                            if (runNumber[arg] == 36)
                            {
                                latex[arg] = new TLatex(peakValleyGraph0_14->GetX()[pointCounter0_14] - 60, peakValleyGraph0_14->GetY()[pointCounter0_14] - 0.5, Form("%d", runPosition[arg]));
                            }
                            else
                            {
                                latex[arg] = new TLatex(peakValleyGraph0_14->GetX()[pointCounter0_14] - 60, peakValleyGraph0_14->GetY()[pointCounter0_14], Form("%d", runPosition[arg]));
                            }
                            peakValleyGraph0_14->GetListOfFunctions()->Add(latex[arg]);
                            latex[arg]->SetTextFont(42);
                            latex[arg]->SetTextSize(.04);
                        }
                        sigmaGraph0_14->SetPoint(pointCounter0_14, nPhotons[arg], stdDevs[arg][8]);
                        sigmaGraph0_14->SetPointError(pointCounter0_14, 0, stdDevErrors[arg][8]);
                        pointCounter0_14 += 1;
                    }
                    else if (runEnergy[arg] == 26)
                    {
                        peakValleyGraph0_26->SetPoint(pointCounter0_26, nPhotons[arg], peakValleyRatio[arg]);
                        peakValleyGraph0_26->SetPointError(pointCounter0_26, 0, 0, TMath::Min(peakValleyRatio[arg] - 1, peakValleyRatioErrors[arg]), peakValleyRatioErrors[arg]);
                        if (nPhotons[arg] >= 1000)
                        {
                            latex[arg] = new TLatex(peakValleyGraph0_26->GetX()[pointCounter0_26] - 60, peakValleyGraph0_26->GetY()[pointCounter0_26], Form("%d", runPosition[arg]));
                            peakValleyGraph0_26->GetListOfFunctions()->Add(latex[arg]);
                            latex[arg]->SetTextFont(42);
                            latex[arg]->SetTextSize(.04);
                        }
                        sigmaGraph0_26->SetPoint(pointCounter0_26, nPhotons[arg], stdDevs[arg][8]);
                        sigmaGraph0_26->SetPointError(pointCounter0_26, 0, stdDevErrors[arg][8]);
                        pointCounter0_26 += 1;
                    }
                    else if (runEnergy[arg] == 52)
                    {
                        peakValleyGraph0_52->SetPoint(pointCounter0_52, nPhotons[arg], peakValleyRatio[arg]);
                        peakValleyGraph0_52->SetPointError(pointCounter0_52, 0, 0, TMath::Min(peakValleyRatio[arg] - 1, peakValleyRatioErrors[arg]), peakValleyRatioErrors[arg]);
                        if (nPhotons[arg] >= 1000)
                        {
                            latex[arg] = new TLatex(peakValleyGraph0_52->GetX()[pointCounter0_52] - 60, peakValleyGraph0_52->GetY()[pointCounter0_52], Form("%d", runPosition[arg]));
                            peakValleyGraph0_52->GetListOfFunctions()->Add(latex[arg]);
                            latex[arg]->SetTextFont(42);
                            latex[arg]->SetTextSize(.04);
                        }
                        sigmaGraph0_52->SetPoint(pointCounter0_52, nPhotons[arg], stdDevs[arg][8]);
                        sigmaGraph0_52->SetPointError(pointCounter0_52, 0, stdDevErrors[arg][8]);
                        pointCounter0_52 += 1;
                    }

                }
            }
            else
            {
                if (beamAngles[arg])
                {
                    if (peakValues[arg][ch] < -50.0 && ch == 8) {
                        peakGraph30[ch]->SetPoint(pointCounter30, beamAngles[arg], peakValues[arg][ch] + 360.0); //Fold up since the fit routine gets the wrong peak sometimes
                    }
                    else {
                        peakGraph30[ch]->SetPoint(pointCounter30, beamAngles[arg], peakValues[arg][ch]);
                    }
                    peakGraph30[ch]->SetPointError(pointCounter30, 0, peakValueErrors[arg][ch]);
                    stdDevGraph30[ch]->SetPoint(pointCounter30, beamAngles[arg], stdDevs[arg][ch]);
                    stdDevGraph30[ch]->SetPointError(pointCounter30, 0, stdDevErrors[arg][ch]);
                    pointCounter30 += 1;
                }
                if (peakValleyRatio[arg])
                {
                    if (runEnergy[arg] == 14)
                    {
                        peakValleyGraph30_14->SetPoint(pointCounter30_14, nPhotons[arg], peakValleyRatio[arg]);
                        peakValleyGraph30_14->SetPointError(pointCounter30_14, 0, 0, TMath::Min(peakValleyRatio[arg] - 1, peakValleyRatioErrors[arg]), peakValleyRatioErrors[arg]);
                        sigmaGraph30_14->SetPoint(pointCounter30_14, nPhotons[arg], stdDevs[arg][8]);
                        sigmaGraph30_14->SetPointError(pointCounter30_14, 0, stdDevErrors[arg][8]);
                        pointCounter30_14 += 1;
                    }
                    else if (runEnergy[arg] == 26)
                    {
                        peakValleyGraph30_26->SetPoint(pointCounter30_26, nPhotons[arg], peakValleyRatio[arg]);
                        peakValleyGraph30_26->SetPointError(pointCounter30_26, 0, 0, TMath::Min(peakValleyRatio[arg] - 1, peakValleyRatioErrors[arg]), peakValleyRatioErrors[arg]);
                        sigmaGraph30_26->SetPoint(pointCounter30_26, nPhotons[arg], stdDevs[arg][8]);
                        sigmaGraph30_26->SetPointError(pointCounter30_26, 0, stdDevErrors[arg][8]);
                        pointCounter30_26 += 1;
                    }
                    else if (runEnergy[arg] == 52)
                    {
                        peakValleyGraph30_52->SetPoint(pointCounter30_52, nPhotons[arg], peakValleyRatio[arg]);
                        peakValleyGraph30_52->SetPointError(pointCounter30_52, 0, 0, TMath::Min(peakValleyRatio[arg] - 1, peakValleyRatioErrors[arg]), peakValleyRatioErrors[arg]);
                        sigmaGraph30_52->SetPoint(pointCounter30_52, nPhotons[arg], stdDevs[arg][8]);
                        sigmaGraph30_52->SetPointError(pointCounter30_52, 0, stdDevErrors[arg][8]);
                        pointCounter30_52 += 1;
                    }

                }
            }            
        }        
    }

    //Filling correlation Graphs:
    int pointCounterCorrelation0 = 0;
    int pointCounterCorrelation30 = 0;

    for (int arg = 1; arg < argc; arg++)
    {
        if (runAngle[arg] == 0)
        {
            correlationGraph0->SetPoint(pointCounterCorrelation0, stdDevs[arg][8], peakValleyRatio[arg]);
            correlationGraph0->SetPointError(pointCounterCorrelation0, stdDevErrors[arg][8], peakValleyRatioErrors[arg]);
            pointCounterCorrelation0 += 1;
        }
        else
        {
            correlationGraph30->SetPoint(pointCounterCorrelation30, stdDevs[arg][8], peakValleyRatio[arg]);
            correlationGraph30->SetPointError(pointCounterCorrelation30, stdDevErrors[arg][8], peakValleyRatioErrors[arg]);
            pointCounterCorrelation30 += 1;
        }
    }
    //FILLING GRAPHS WITH POINTS DONE
    
    //Setting up graph style:
    peakValleyGraph0_14->SetMarkerColor(kBlue - 7);
    peakValleyGraph0_14->SetMarkerStyle(47);
    peakValleyGraph0_14->SetMarkerSize(2);
    sigmaGraph0_14->SetMarkerColor(kBlue - 7);
    sigmaGraph0_14->SetMarkerStyle(47);
    sigmaGraph0_14->SetMarkerSize(2);
    //peakValleyFit0_14->SetLineWidth(3);
    //peakValleyFit0_14->SetLineColor(kGreen + 3);
    peakValleyGraph0_26->SetMarkerColor(kGreen + 2);
    peakValleyGraph0_26->SetMarkerStyle(47);
    peakValleyGraph0_26->SetMarkerSize(2);
    sigmaGraph0_26->SetMarkerColor(kGreen + 2);
    sigmaGraph0_26->SetMarkerStyle(47);
    sigmaGraph0_26->SetMarkerSize(2);
    //peakValleyFit0_26->SetLineWidth(3);
    //peakValleyFit0_26->SetLineColor(kGreen + 2);
    peakValleyGraph0_52->SetMarkerColor(kRed - 6);
    peakValleyGraph0_52->SetMarkerStyle(47);
    peakValleyGraph0_52->SetMarkerSize(2);
    sigmaGraph0_52->SetMarkerColor(kRed - 6);
    sigmaGraph0_52->SetMarkerStyle(47);
    sigmaGraph0_52->SetMarkerSize(2);
    //peakValleyFit0_52->SetLineWidth(3);
    //peakValleyFit0_52->SetLineColor(kGreen + 1);

    peakValleyGraph30_14->SetMarkerColor(kBlue - 7);
    peakValleyGraph30_14->SetMarkerStyle(47);
    peakValleyGraph30_14->SetMarkerSize(2);
    sigmaGraph30_14->SetMarkerColor(kBlue - 7);
    sigmaGraph30_14->SetMarkerStyle(47);
    sigmaGraph30_14->SetMarkerSize(2);
    //peakValleyFit30_14->SetLineWidth(3);
    //peakValleyFit30_14->SetLineColor(kRed - 1);
    peakValleyGraph30_26->SetMarkerColor(kGreen + 2);
    peakValleyGraph30_26->SetMarkerStyle(47);
    peakValleyGraph30_26->SetMarkerSize(2);
    sigmaGraph30_26->SetMarkerColor(kGreen + 2);
    sigmaGraph30_26->SetMarkerStyle(47);
    sigmaGraph30_26->SetMarkerSize(2);
    //peakValleyFit30_26->SetLineWidth(3);
    //peakValleyFit30_26->SetLineColor(kRed - 2);
    peakValleyGraph30_52->SetMarkerColor(kRed - 6);
    peakValleyGraph30_52->SetMarkerStyle(47);
    peakValleyGraph30_52->SetMarkerSize(2);
    sigmaGraph30_52->SetMarkerColor(kRed - 6);
    sigmaGraph30_52->SetMarkerStyle(47);
    sigmaGraph30_52->SetMarkerSize(2);
    //peakValleyFit30_52->SetLineWidth(3);
    //peakValleyFit30_52->SetLineColor(kRed - 3);

    correlationGraph0->SetMarkerColor(kGreen + 1);
    correlationGraph0->SetMarkerStyle(20);
    correlationGraph0->SetMarkerSize(2);
    correlationGraph30->SetMarkerColor(kBlue - 2);
    correlationGraph30->SetMarkerStyle(20);
    correlationGraph30->SetMarkerSize(2);

    

    for (int i = 0; i < nCh + 1; i++)
    {
        peakGraph0[i]->SetMarkerStyle(20);
        peakGraph0[i]->SetMarkerSize(1.7);
        peakGraph0[i]->SetMarkerColor(kGreen + 1);
        peakGraph30[i]->SetMarkerStyle(34);
        peakGraph30[i]->SetMarkerSize(1.7);
        peakGraph30[i]->SetMarkerColor(kBlue - 2);
        stdDevGraph0[i]->SetMarkerStyle(20);
        stdDevGraph0[i]->SetMarkerSize(1.7);
        stdDevGraph0[i]->SetMarkerColor(kGreen + 1);
        stdDevGraph30[i]->SetMarkerStyle(34);
        stdDevGraph30[i]->SetMarkerSize(1.7);
        stdDevGraph30[i]->SetMarkerColor(kBlue - 2);
        if (i == 8)
        {
            peakGraph0[i]->SetMarkerSize(3);
            peakGraph30[i]->SetMarkerSize(3);
            stdDevGraph0[i]->SetMarkerSize(3);
            stdDevGraph30[i]->SetMarkerSize(3);
        }
    }

    //Creating legends:
    for (int i = 0; i < nCh + 1; i++)
    {
        peakLegend[i] = new TLegend(peakLegendXMin[i], peakLegendYMin[i], peakLegendXMax[i], peakLegendYMax[i]);
        peakLegend[i]->SetFillColorAlpha(kWhite, 0.3);
        if (i != 8)
        {
            stdDevLegend[i] = new TLegend(0.55, 0.65, 0.935, 0.9);
        }
        else
        {
            stdDevLegend[i] = new TLegend(0.54, 0.65, 0.9, 0.9);
        }
        
        stdDevLegend[i]->SetFillColorAlpha(kWhite, 0.7);
        peakLegend[i]->AddEntry(peakGraph0[i], "0 Degrees", "PE");
        peakLegend[i]->AddEntry(peakGraph30[i], "30 Degrees", "PE");
        stdDevLegend[i]->AddEntry(stdDevGraph0[i], "0 Degrees", "PE");
        stdDevLegend[i]->AddEntry(stdDevGraph30[i], "30 Degrees", "PE");
    }
    peakValleyLegend0->AddEntry(peakValleyGraph0_14, "0 Degrees, 1.4 GeV", "PE");
    peakValleyLegend0->AddEntry(peakValleyGraph0_26, "0 Degrees, 2.6 GeV", "PE");
    peakValleyLegend0->AddEntry(peakValleyGraph0_52, "0 Degrees, 5.2 GeV", "PE");
    peakValleyLegend30->AddEntry(peakValleyGraph30_14, "30 Degrees, 1.4 GeV", "PE");
    peakValleyLegend30->AddEntry(peakValleyGraph30_26, "30 Degrees, 2.6 GeV", "PE");
    peakValleyLegend30->AddEntry(peakValleyGraph30_52, "30 Degrees, 5.2 GeV", "PE");
    peakValleyLegend0->SetFillColorAlpha(kWhite, 0.7);
    peakValleyLegend30->SetFillColorAlpha(kWhite, 0.7);

    sigmaLegend0->AddEntry(sigmaGraph0_14, "0 Degrees, 1.4 GeV", "PE");
    sigmaLegend0->AddEntry(sigmaGraph0_26, "0 Degrees, 2.6 GeV", "PE");
    sigmaLegend0->AddEntry(sigmaGraph0_52, "0 Degrees, 5.2 GeV", "PE");
    sigmaLegend30->AddEntry(sigmaGraph30_14, "30 Degrees, 1.4 GeV", "PE");
    sigmaLegend30->AddEntry(sigmaGraph30_26, "30 Degrees, 2.6 GeV", "PE");
    sigmaLegend30->AddEntry(sigmaGraph30_52, "30 Degrees, 5.2 GeV", "PE");
    sigmaLegend0->SetFillColorAlpha(kWhite, 0.7);
    sigmaLegend30->SetFillColorAlpha(kWhite, 0.7);

    correlationLegend->AddEntry(correlationGraph0, "0 Degrees", "PE");
    correlationLegend->AddEntry(correlationGraph30, "30 Degrees", "PE");
    correlationLegend->SetFillColorAlpha(kWhite, 0.7);

    //Merging and printing graphs with legends:
    for (int i = 0; i < nCh + 1; i++)
    {
        if (i == 8)
        {
            peakMultiGraph[i] = new TMultiGraph("Peak Positions for all Channels Included", "Peak Positions for all Channels Included;#phi_{beam} [deg.];#phi_{peak} [deg.]");
            stdDevMultiGraph[i] = new TMultiGraph("Standard Deviations for all Channels Included", "Standard Deviations for all Channels Included;#phi_{beam} [deg.];Central Peak #sigma");
        }
        else
        {
            peakMultiGraph[i] = new TMultiGraph(Form("Ignoring Channel %d (%d deg.);#phi_{beam} [deg.];#phi_{peak} [deg.]", i, angleList[i]),
                Form("Ignoring Channel %d (%d deg.)", i, angleList[i])); //;#phi_{beam} [deg.];#phi_{peak} [deg.]
            stdDevMultiGraph[i] = new TMultiGraph(Form("Ignoring Channel %d (%d deg.);#phi_{beam} [deg.];#sigma", i, angleList[i]),
                Form("Ignoring Channel %d (%d deg.)", i, angleList[i])); //;#phi_{beam} [deg.];#sigma
        }

        peakMultiGraph[i]->Add(peakGraph0[i]);
        peakMultiGraph[i]->Add(peakGraph30[i]);
        peakMultiGraph[i]->GetXaxis()->SetTitleSize(.07);
        peakMultiGraph[i]->GetXaxis()->SetLabelSize(.06);
        peakMultiGraph[i]->GetYaxis()->SetTitleSize(.07);
        peakMultiGraph[i]->GetYaxis()->SetLabelSize(.06);

        stdDevMultiGraph[i]->Add(stdDevGraph0[i]);
        stdDevMultiGraph[i]->Add(stdDevGraph30[i]);
        stdDevMultiGraph[i]->GetXaxis()->SetTitleSize(.07);
        stdDevMultiGraph[i]->GetXaxis()->SetLabelSize(.06);
        stdDevMultiGraph[i]->GetYaxis()->SetTitleSize(.07);
        stdDevMultiGraph[i]->GetYaxis()->SetLabelSize(.06);

        if (i != 8)
        {
            peakPad.cd(i + 1);
            gPad->SetGrid(20, 2);
            gPad->SetLeftMargin(.065);
            gPad->SetBottomMargin(.052);
            gPad->SetRightMargin(.065);
            peakMultiGraph[i]->Draw("APE");
            peakLegend[i]->Draw("same");
            stdDevPad.cd(i + 1);
            gPad->SetGrid(20, 10);
            gPad->SetLeftMargin(.065);
            gPad->SetBottomMargin(.052);
            gPad->SetRightMargin(.065);
            stdDevMultiGraph[i]->Draw("APE");
            stdDevLegend[i]->Draw("same");
        }
        else
        {
            peakCanvas8.cd();
            gPad->SetGrid(20, 2);
            gPad->SetLeftMargin(.18);
            gPad->SetBottomMargin(.16);
            peakMultiGraph[i]->Draw("APE");
            peakLegend[i]->Draw("same");
            stdDevCanvas8.cd();
            gPad->SetGrid(20, 10);
            gPad->SetLeftMargin(.18);
            gPad->SetBottomMargin(.16);
            stdDevMultiGraph[i]->Draw("APE");
            stdDevLegend[i]->Draw("same");
        }
        
    }

    peakValleyMultiGraph0->Add(peakValleyGraph0_14);
    peakValleyMultiGraph0->Add(peakValleyGraph0_26);
    peakValleyMultiGraph0->Add(peakValleyGraph0_52);
    peakValleyMultiGraph30->Add(peakValleyGraph30_14);
    peakValleyMultiGraph30->Add(peakValleyGraph30_26);
    peakValleyMultiGraph30->Add(peakValleyGraph30_52);

    sigmaMultiGraph0->Add(sigmaGraph0_14);
    sigmaMultiGraph0->Add(sigmaGraph0_26);
    sigmaMultiGraph0->Add(sigmaGraph0_52);
    sigmaMultiGraph30->Add(sigmaGraph30_14);
    sigmaMultiGraph30->Add(sigmaGraph30_26);
    sigmaMultiGraph30->Add(sigmaGraph30_52);

    correlationMultiGraph->Add(correlationGraph0);
    correlationMultiGraph->Add(correlationGraph30);

    peakValleyCanvas0.cd(0);
    gPad->SetLeftMargin(.18);
    gPad->SetBottomMargin(.16);
    peakValleyCanvas0.SetGrid();
    peakValleyMultiGraph0->SetTitle("Maximum to Minimum Ratio #it{r} for all Channels for 0 deg.;Average light yield per Event;#it{r}");
    peakValleyMultiGraph0->GetXaxis()->SetTitleSize(.07);
    peakValleyMultiGraph0->GetXaxis()->SetLabelSize(.06);
    peakValleyMultiGraph0->GetYaxis()->SetTitleSize(.07);
    peakValleyMultiGraph0->GetYaxis()->SetLabelSize(.06);
    peakValleyMultiGraph0->Draw("AP");
    peakValleyLegend0->Draw("same");

    sigmaCanvas0.cd(0);
    gPad->SetLeftMargin(.18);
    gPad->SetBottomMargin(.16);
    sigmaCanvas0.SetGrid(200, 20);
    sigmaMultiGraph0->SetTitle("Comparison of #sigma values for all Channels for 0 deg.;Average light yield per Event;#sigma");
    sigmaMultiGraph0->Draw("AP");
    sigmaMultiGraph0->GetXaxis()->SetTitleSize(.07);
    sigmaMultiGraph0->GetXaxis()->SetLabelSize(.06);
    sigmaMultiGraph0->GetYaxis()->SetTitleSize(.07);
    sigmaMultiGraph0->GetYaxis()->SetLabelSize(.06);
    sigmaLegend0->Draw("same");

    peakValleyCanvas30.cd(0);
    gPad->SetLeftMargin(.18);
    gPad->SetBottomMargin(.16);
    peakValleyCanvas30.SetGrid();
    peakValleyMultiGraph30->SetTitle("Maximum to Minimum Ratio #it{r} for All Channels for 30 deg.;Average Number of Photoelectrons per Event;#it{r}");
    peakValleyMultiGraph30->GetXaxis()->SetTitleSize(.07);
    peakValleyMultiGraph30->GetXaxis()->SetLabelSize(.06);
    peakValleyMultiGraph30->GetYaxis()->SetTitleSize(.07);
    peakValleyMultiGraph30->GetYaxis()->SetLabelSize(.06);
    peakValleyMultiGraph30->Draw("AP");
    peakValleyLegend30->Draw("same");

    sigmaCanvas30.cd(0);
    gPad->SetLeftMargin(.18);
    gPad->SetBottomMargin(.16);
    sigmaCanvas30.SetGrid(100, 50);
    sigmaMultiGraph30->SetTitle("Comparison of #sigma values for all Channels for 30 deg.;Average Number of Photoelectrons per Event;#sigma");
    sigmaMultiGraph30->Draw("AP");
    sigmaMultiGraph30->GetXaxis()->SetTitleSize(.07);
    sigmaMultiGraph30->GetXaxis()->SetLabelSize(.06);
    sigmaMultiGraph30->GetYaxis()->SetTitleSize(.07);
    sigmaMultiGraph30->GetYaxis()->SetLabelSize(.06);
    sigmaLegend30->Draw("same");

    canvas91011.cd();
    gPad->SetGrid();
    legend91011->Draw("same");

    correlationCanvas.cd(0);
    gPad->SetGrid();
    gPad->SetLeftMargin(.18);
    gPad->SetBottomMargin(.16);
    correlationMultiGraph->GetXaxis()->SetTitleSize(.07);
    correlationMultiGraph->GetXaxis()->SetLabelSize(.06);
    correlationMultiGraph->GetYaxis()->SetTitleSize(.07);
    correlationMultiGraph->GetYaxis()->SetLabelSize(.06);
    correlationMultiGraph->Draw("APE");
    correlationLegend->Draw("same");

    peakCanvas.SaveAs(Form("%speakValues.pdf", outPath.c_str()));
    peakCanvas8.SaveAs(Form("%speakValuesAllCh.pdf", outPath.c_str()));
    stdDevCanvas.SaveAs(Form("%sstdDevs.pdf", outPath.c_str()));
    stdDevCanvas8.SaveAs(Form("%sstdDevsAllCh.pdf", outPath.c_str()));
    correlationCanvas.SaveAs(Form("%sCorrelations.pdf", outPath.c_str()));

    peakValleyCanvas0.SaveAs(Form("%speakValleyRatios0deg.pdf", outPath.c_str()));
    peakValleyCanvas30.SaveAs(Form("%speakValleyRatios30deg.pdf", outPath.c_str()));
    sigmaCanvas0.SaveAs(Form("%ssigmas0deg.pdf", outPath.c_str()));
    sigmaCanvas30.SaveAs(Form("%ssigmas30deg.pdf", outPath.c_str()));

    canvas91011.SaveAs(Form("%spositionComparison.pdf", outPath.c_str()));
}