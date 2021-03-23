#include <TGraph.h>
#include <TGraphErrors.h>
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
#include <TEllipse.h>
#include <TPaveLabel.h>
#include <TLatex.h>

//Including standard C/C++-Libraries:
#include <sys/stat.h>
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

using namespace TMath;


//###############################################
//########## DEFINING THE PHOTON CLASS ##########
//###############################################

class Photon
{
private:    //Attributes:
    float womX = -0.07;
    float womY = 0.07;
    float womRadius = 0.03;     //Central point and radius of WOM D
    float angle;                //Angle under which the scintillation photon is emitted from the beam position (where 0 is positive x direction). Not to be confused with phi
    float xPos, yPos;           //Coordinates of beam positions (all numbers in meters) in the cartesian coordinate system with the ORIGIN AT THE BOX CORNER
    float stepLength = 0.0001;
    float distanceWalked = 0;
    float phi;               //Coordinate of the scintillation photon in the polar coordinate system with the ORIGIN AT THE CENTER OF WOM D. Here, phi=0 is negative x direction
    float attLength = 0.0853; //attentuation length from Patricks measurements in m (for 360nm light)

public:     //Member functions:
    void setBeamPosition(int pos)
    {
        switch (pos)
        {
            case 2:
                xPos = -0.06;
                yPos = 0.90;
                break;
            case 9:
                xPos = -0.22;
                yPos = 0.07;
                break;
            case 10:
                xPos = -0.176;
                yPos = 0.176;
                break;
            case 11:
                xPos = -0.07;
                yPos = 0.22;
                break;
        }
    }
    void setAngle(float newAngle) //Sets the emission angle at the beam position
    {
        angle = newAngle;
    }
    float getAngle()
    {
        return angle;
    }
    void incrementStep()        //Photon walks a step of 1mm length in the direction specified by angle
    {
        xPos += stepLength * Cos(angle);
        yPos += stepLength * Sin(angle);
        distanceWalked += stepLength;
    }
    bool isInWom()              //Checks whether the photon has entered the wom
    {
        if (Sqrt(Power(xPos - womX, 2) + Power(yPos - womY, 2)) <= womRadius)
        {
            return true;
        }
        else {
            return false;
        }
    }
    void reflectIfOutsideBox()  //Changes the photon step direction if the photon leaves the box
    {
        if (xPos > 0)
        {
            angle = Pi() - angle;
        }
        if (yPos < 0)
        {
            angle = -angle;
        }
    }
    float getXPos()
    {
        return xPos;
    }
    float getYPos()
    {
        return yPos;
    }
    int getPhi()
    {
        if (xPos < womX)
        {
            phi = Floor(((ATan((yPos - 0.07) / (xPos + 0.07))) - Pi()) * 180.0 / Pi());
        }
        else
        {
            phi = Floor((ATan((yPos - 0.07) / (xPos + 0.07))) * 180.0 / Pi());   //Returns angle in Degrees, rounded down
        }
        if (phi >= 0)
        {
            return phi;
        }
        else
        {
            return phi + 360;
        }
    }
    float getIntensity()
    {
        return (1.0 / distanceWalked) * Exp(-distanceWalked / attLength);
    }
};

//###################################
//########## MAIN FUNCTION ##########
//###################################

int main()
{
    //Simulation is done twice, once with reflection and once without. Number 2 at the end of objects <--> no reflection!

    std::string out_path = "../runs/WomDSimulation/";
    const int numberOfPositions = 3;                        //Currently positions 9,10,11
    int positionVec[] = { 9, 10, 11 };                   //In the way the code is currently written, the reference position needs to be at positionVec[1]!
    const int numberOfDataPoints = 8;                       //Corresponds to 8 channels, but can in principle be chosen arbitrarily
    const int numberOfAngleDivisions = 40;                  //Precision; a photon is sent out every 1 / numberOfAngleDivisions degrees
    const float fNumberOfAngleDivisions = numberOfAngleDivisions * 1.0;//Same value but as a floating point number, to be used everywhere where calculations with float results need to be made
    const int numberOfAngles = 360 * numberOfAngleDivisions;//Possible number of angles under which a photon can be sent out.
    const int numberOfSteps = 10000;                        //Maximum number of steps a photon takes
    const int colorList[numberOfPositions] = { 433,617,797 };           //Color list out of which the colors for drawing the simulation and the intensities are taken: Orange, Red, Green, Blue
                                                            //Check https://root.cern.ch/doc/master/classTColor.html for details
    const int markerList[numberOfPositions] = { 20,21,34 };              //Marker List for the same purpose; https://root.cern.ch/doc/master/pict1_TAttMarker_002.png
    int sectionWidth = 360 / numberOfDataPoints;
    int halfSectionWidth = 360 / (2 * numberOfDataPoints);
    gStyle->SetLineScalePS(1);                              //Resolution
    gStyle->SetHatchesSpacing(0.5);                         //Spacing between lines in boxes in legends
    gStyle->SetGridStyle(3);
    gStyle->SetGridWidth(1);
    gStyle->SetGridColor(16);
    gStyle->SetTitleSize(0.08, "t");

    //Canvases:
    TCanvas* c1 = new TCanvas("c1", "WOM D Simulation", 1080, 2700);
    TCanvas* c2 = new TCanvas("c2", "Simulated Intensities", 1920, 1080);
    TCanvas* c3 = new TCanvas("c3", "Simulated Reduced Intensities", 1920, 1080);
    TCanvas* c4 = new TCanvas("c4", "Simulation without Reflection", 1080, 2700);
    TCanvas* c5 = new TCanvas("c5", "Simulated Intensities without Reflection", 1920, 1080);
    TCanvas* c6 = new TCanvas("c6", "Simulated #it{I}_{red} without Reflection", 1920, 1080);

    c1->cd();
    TPaveLabel title1(0.1, 0.96, 0.9, 0.99, "WOM D Simulation");
    title1.SetTextSize(.55);
    title1.SetLineColor(0);
    title1.SetBorderSize(0);
    title1.SetFillColor(0);
    title1.Draw();
    TPad graphPad1("Graphs", "Graphs", 0, 0, 1, 0.95);
    graphPad1.Draw();
    graphPad1.cd();
    graphPad1.Divide(1, 3, 0.001, 0.01);

    c4->cd();
    TPaveLabel title4(0.1, 0.96, 0.9, 0.99, "WOM D Simulation without Reflection");
    title4.SetTextSize(.55);
    title4.SetLineColor(0);
    title4.SetBorderSize(0);
    title4.SetFillColor(0);
    title4.Draw();
    TPad graphPad4("Graphs", "Graphs", 0, 0, 1, 0.95);
    graphPad4.Draw();
    graphPad4.cd();
    graphPad4.Divide(1, 3, 0.001, 0.01);

    //Graphs for drawing of simulation with reflection:
    TGraph* Gr[numberOfPositions][numberOfAngles / numberOfAngleDivisions];     //Increasing first index corresponds to increasing POSITION number (here 0-> 2, 1->9, 2->10, 3->11)
    TGraph* GrDummy[numberOfPositions];                      //Add points to be included in the multigraph, so that the axis range is set right; SetLimits() does not work for some reason
    TMultiGraph* PosGraph[numberOfPositions];
    TLegend* PosLegend[numberOfPositions];

    //And without reflection:
    TGraph* Gr2[numberOfPositions][numberOfAngles / numberOfAngleDivisions];
    TGraph* GrDummy2[numberOfPositions];
    TMultiGraph* PosGraph2[numberOfPositions];
    TLegend* PosLegend2[numberOfPositions];

    //Other Objects for drawing of simulation:
    TEllipse* womCircle = new TEllipse(-0.07, 0.07, 0.03, 0.03);
    womCircle->SetLineWidth(5);
    //TLatex* womLabel10 = new TLatex(-0.09, 0.065, "WOM D");
    TLatex* womLabel = new TLatex(-0.0895, 0.064, "WOM D");

    //Graphs for intensity results with reflection:
    TGraph* PosIntGraph[numberOfPositions];
    TGraph* PosIntGraphNorm[numberOfPositions];
    TMultiGraph* IntGraph = new TMultiGraph();
    TMultiGraph* IntGraphNorm = new TMultiGraph();
    TLegend* IntLegend = new TLegend(0.65, 0.7, 0.9, 0.9);
    IntLegend->SetFillColorAlpha(kWhite, .88);
    TLegend* IntLegendNorm = new TLegend(0.65, 0.7, 0.9, 0.9);
    IntLegendNorm->SetFillColorAlpha(kWhite, .8);

    //And without reflection:
    TGraph* PosIntGraph2[numberOfPositions];
    TGraph* PosIntGraphNorm2[numberOfPositions];
    TMultiGraph* IntGraph2 = new TMultiGraph();
    TMultiGraph* IntGraphNorm2 = new TMultiGraph();
    TLegend* IntLegend2 = new TLegend(0.65, 0.7, 0.9, 0.9);
    IntLegend2->SetFillColorAlpha(kWhite, .88);
    TLegend* IntLegendNorm2 = new TLegend(0.65, 0.7, 0.9, 0.9);
    IntLegendNorm2->SetFillColorAlpha(kWhite, .8);

    //Variables containing intensity results:
    float intensities[numberOfPositions][numberOfDataPoints];
    float intensitiesNorm[numberOfPositions][numberOfDataPoints];
    float totalInt[numberOfPositions];
    float intensities2[numberOfPositions][numberOfDataPoints];
    float intensitiesNorm2[numberOfPositions][numberOfDataPoints];
    float totalInt2[numberOfPositions];

    //Variable for setting corresponding graph points in simulation (is counted up during one loop and reset at the beginning of the next):
    int pointNumber;

    //Photon object arrays:
    Photon Ph[numberOfPositions][numberOfAngles];
    Photon Ph2[numberOfPositions][numberOfAngles];

    //With reflection:
    //Looping over all positions:
    for (int i = 0; i < numberOfPositions; i++)
    {
        PosGraph[i] = new TMultiGraph();

        //Looping over all angles:
        for (int k = 0; k < numberOfAngles; k++)
        {
            //For each photon the starting position and the starting angle is set:
            Ph[i][k].setBeamPosition(positionVec[i]);
            Ph[i][k].setAngle((k / fNumberOfAngleDivisions) * Pi() / 180.0);
            pointNumber = 0;

            //The photon trajectory is calculated for every angle, but a graph is drawn only once every degree, to improve optical representation and reduce file size:
            if (k % numberOfAngleDivisions == 0)
            {
                Gr[i][k / numberOfAngleDivisions] = new TGraph();
                for (int step = 0; step < numberOfSteps; step++)
                {
                    if (step % 10 == 0)
                    {
                        //For every 10th step, a point is added to the graph:
                        Gr[i][k / numberOfAngleDivisions]->SetPoint(pointNumber, Ph[i][k].getXPos(), Ph[i][k].getYPos());
                        pointNumber++;
                    }
                    //Check whether the photon hit the wom or a wall, and in the end, the photon walks another step :
                    if (Ph[i][k].isInWom())
                    {
                        break;
                    }
                    Ph[i][k].reflectIfOutsideBox();
                    Ph[i][k].incrementStep();
                }
                if (Gr[i][k / numberOfAngleDivisions]->GetN() < numberOfSteps / 10)
                {
                    //Only draw the graph if it has less than numberOfSteps points (meaning the photon hit the WOM before finishing the loop). Otherwise, delete it:
                    Gr[i][k / numberOfAngleDivisions]->SetMarkerColorAlpha(colorList[i] - 3, .8);
                    Gr[i][k / numberOfAngleDivisions]->SetLineColor(colorList[i] - 3);
                    Gr[i][k / numberOfAngleDivisions]->SetFillColor(colorList[i] - 3);
                    PosGraph[i]->Add(Gr[i][k / numberOfAngleDivisions]);
                }
                else
                {
                    Gr[i][k / numberOfAngleDivisions]->Delete();
                }
            }
            else
            {
                for (int step = 0; step < numberOfSteps; step++)
                {
                    if (Ph[i][k].isInWom())
                    {
                        break;
                    }
                    Ph[i][k].reflectIfOutsideBox();
                    Ph[i][k].incrementStep();
                }
            }
            int deg = Ph[i][k].getPhi();                     //This retrieves the rounded down angle at which the scintillation photon hits the WOM tube       
            int section;                                     //The section in which the photon hits the WOM

            //Dividing the WOM Tube into numberOfDataPoints sections, such that 0 degrees is in the middle of a section, and checking whether the photon hit the WOM in each section:
            if ((0 <= deg && deg < halfSectionWidth) || ((360 - halfSectionWidth) <= deg && deg < 360))
            {
                section = 0;
            }
            else
            {
                for (int div = 1; div < numberOfDataPoints; div++)
                {
                    if (halfSectionWidth + (div - 1) * sectionWidth <= deg && deg < halfSectionWidth + div * sectionWidth)
                    {
                        section = div;
                    }
                }
            }
            intensities[i][section] += Ph[i][k].getIntensity();       //This takes the intensity at that angle and adds it to the entry belonging to the angle in the intensities array
        }   //End of loop over all angles

        //Calculating the total intensity for each position and dividing each data point by it:
        for (int div = 0; div < numberOfDataPoints; div++)
        {
            totalInt[i] += intensities[i][div];
        }
        for (int div = 0; div < numberOfDataPoints; div++)
        {
            intensities[i][div] = intensities[i][div] / totalInt[i];
        }
        printf("Simulation of Position %d finished. Total intensity:  \t%f\n", positionVec[i], totalInt[i]);
    }   //End of loop over all positions

    //Without reflection:
    //Looping over all positions:
    for (int i = 0; i < numberOfPositions; i++)
    {
        PosGraph2[i] = new TMultiGraph();

        //Looping over all angles:
        for (int k = 0; k < numberOfAngles; k++)
        {
            //For each photon the starting position and the starting angle is set:
            Ph2[i][k].setBeamPosition(positionVec[i]);
            Ph2[i][k].setAngle((k / fNumberOfAngleDivisions) * Pi() / 180.0);
            pointNumber = 0;

            //The photon trajectory is calculated for every angle, but a graph is drawn only once every degree, to improve optical representation and reduce file size:
            if (k % numberOfAngleDivisions == 0)
            {
                Gr2[i][k / numberOfAngleDivisions] = new TGraph();
                for (int step = 0; step < numberOfSteps; step++)
                {
                    if (step % 10 == 0)
                    {
                        //For every 10th step, a point is added to the graph:
                        Gr2[i][k / numberOfAngleDivisions]->SetPoint(pointNumber, Ph2[i][k].getXPos(), Ph2[i][k].getYPos());
                        pointNumber++;
                    }
                    if (Ph2[i][k].isInWom())
                    {
                        break;
                    }
                    //NO REFLECTION HERE!
                    Ph2[i][k].incrementStep();
                }
                if (Gr2[i][k / numberOfAngleDivisions]->GetN() < numberOfSteps / 10)
                {
                    //Only draw the graph if it has less than numberOfSteps points (meaning the photon hit the WOM before finishing the loop). Otherwise, delete it:
                    Gr2[i][k / numberOfAngleDivisions]->SetMarkerColorAlpha(colorList[i] - 3, .8);
                    Gr2[i][k / numberOfAngleDivisions]->SetLineColor(colorList[i] - 3);
                    Gr2[i][k / numberOfAngleDivisions]->SetFillColor(colorList[i] - 3);
                    PosGraph2[i]->Add(Gr2[i][k / numberOfAngleDivisions]);
                }
                else
                {
                    Gr2[i][k / numberOfAngleDivisions]->Delete();
                }
            }
            else
            {
                for (int step = 0; step < numberOfSteps; step++)
                {
                    if (Ph2[i][k].isInWom())
                    {
                        break;
                    }
                    //NO REFLECTION HERE!
                    Ph2[i][k].incrementStep();
                }
            }
            int deg = Ph2[i][k].getPhi();                    //This retrieves the rounded down angle at which the scintillation photon hits the WOM tube       
            int section;                                     //The section in which the photon hits the WOM

            //Dividing the WOM Tube into numberOfDataPoints sections, such that 0 degrees is in the middle of a section, and checking whether the photon hit the WOM in each section:
            if ((0 <= deg && deg < halfSectionWidth) || ((360 - halfSectionWidth) <= deg && deg < 360))
            {
                section = 0;
            }
            else
            {
                for (int div = 1; div < numberOfDataPoints; div++)
                {
                    if (halfSectionWidth + (div - 1) * sectionWidth <= deg && deg < halfSectionWidth + div * sectionWidth)
                    {
                        section = div;
                    }
                }
            }
            intensities2[i][section] += Ph2[i][k].getIntensity();       //This takes the intensity at that angle and adds it to the entry belonging to the angle in the intensities array
        }   //End of loop over all angles

        //Calculating the total intensity for each position and dividing each data point by it:
        for (int div = 0; div < numberOfDataPoints; div++)
        {
            totalInt2[i] += intensities2[i][div];
        }
        for (int div = 0; div < numberOfDataPoints; div++)
        {
            intensities2[i][div] = intensities2[i][div] / totalInt2[i];
        }
        printf("Simulation of Position %d without Reflection finished. Total intensity:  \t%f\n", positionVec[i], totalInt2[i]);
    }   //End of loop over all positions

    //Dividing all data points by values for position 10:
    for (int j = 0; j < numberOfDataPoints; j++)
    {
        intensitiesNorm[1][j] = 1;
        intensitiesNorm2[1][j] = 1;
        if (Abs(intensities[1][j]) < Power(10, -3))  //Don't divide by zero! 10^(-6) was arbitrarily chosen but should be fine. If there is actually an entry, it won't be that small
        {
            intensitiesNorm[0][j] = 0;
            intensitiesNorm[2][j] = 0;
        }
        else
        {
            intensitiesNorm[0][j] = intensities[0][j] / intensities[1][j];
            intensitiesNorm[2][j] = intensities[2][j] / intensities[1][j];
        }
        if (Abs(intensities2[1][j]) < Power(10, -3)) 
        {
            intensitiesNorm2[0][j] = 0;
            intensitiesNorm2[2][j] = 0;
        }
        else
        {
            intensitiesNorm2[0][j] = intensities2[0][j] / intensities2[1][j];
            intensitiesNorm2[2][j] = intensities2[2][j] / intensities2[1][j];
        }
    }

    //Creating graphs for intensities and divided intensities:
    for (int i = 0; i < numberOfPositions; i++)
    {
        PosIntGraph[i] = new TGraph();
        PosIntGraph2[i] = new TGraph();
        PosIntGraphNorm[i] = new TGraph();
        PosIntGraphNorm2[i] = new TGraph();
        for (int j = 0; j < numberOfDataPoints; j++)
        {
            PosIntGraph[i]->SetPoint(j, j * sectionWidth, intensities[i][j]);
            PosIntGraphNorm[i]->SetPoint(j, j * sectionWidth, intensitiesNorm[i][j]);
            PosIntGraph[i]->SetLineColor(colorList[i] - 3);
            PosIntGraph[i]->SetLineWidth(5);
            PosIntGraph[i]->SetMarkerSize(3);
            PosIntGraph[i]->SetMarkerStyle(47);
            PosIntGraphNorm[i]->SetLineColor(colorList[i] - 3);
            PosIntGraphNorm[i]->SetLineWidth(5);
            PosIntGraphNorm[i]->SetMarkerSize(3);
            PosIntGraphNorm[i]->SetMarkerStyle(47);

            PosIntGraph2[i]->SetPoint(j, j* sectionWidth, intensities2[i][j]);
            PosIntGraphNorm2[i]->SetPoint(j, j* sectionWidth, intensitiesNorm2[i][j]);
            PosIntGraph2[i]->SetLineColor(colorList[i] - 3);
            PosIntGraph2[i]->SetLineWidth(5);
            PosIntGraph2[i]->SetMarkerSize(3);
            PosIntGraph2[i]->SetMarkerStyle(47);
            PosIntGraphNorm2[i]->SetLineColor(colorList[i] - 3);
            PosIntGraphNorm2[i]->SetLineWidth(5);
            PosIntGraphNorm2[i]->SetMarkerSize(3);
            PosIntGraphNorm2[i]->SetMarkerStyle(47);
        }
        IntGraph->Add(PosIntGraph[i]);
        IntLegend->AddEntry(PosIntGraph[i], Form("Position %d", positionVec[i]));
        IntGraphNorm->Add(PosIntGraphNorm[i]);
        IntLegendNorm->AddEntry(PosIntGraphNorm[i], Form("Position %d", positionVec[i]));

        IntGraph2->Add(PosIntGraph2[i]);
        IntLegend2->AddEntry(PosIntGraph2[i], Form("Position %d", positionVec[i]));
        IntGraphNorm2->Add(PosIntGraphNorm2[i]);
        IntLegendNorm2->AddEntry(PosIntGraphNorm2[i], Form("Position %d", positionVec[i]));
    }
    IntGraph->SetTitle(";Azimutal Angle #phi [deg.];Intensity #it{I} [arb. Units]");
    IntGraph->GetXaxis()->SetLimits(0, 360);
    IntGraph->GetXaxis()->SetNdivisions(8, 8, 4, false);        //Setting the x Axis ticks such that there is a tick every 45 degrees
    IntGraph->GetYaxis()->SetRangeUser(0, 0.5);

    IntGraphNorm->SetTitle(";Azimutal Angle #phi [deg.];Reduced Intensity #it{I}_{red}");
    IntGraphNorm->GetXaxis()->SetLimits(0, 360);
    IntGraphNorm->GetXaxis()->SetNdivisions(8, 8, 4, false);
    IntGraphNorm->GetYaxis()->SetRangeUser(0, 10);

    IntGraph2->SetTitle(";Azimutal Angle #phi [deg.];Intensity #it{I} [arb. Units]");
    IntGraph2->GetXaxis()->SetLimits(0, 360);
    IntGraph2->GetXaxis()->SetNdivisions(8, 8, 4, false);        //Setting the x Axis ticks such that there is a tick every 45 degrees
    IntGraph2->GetYaxis()->SetRangeUser(0, 0.55);

    IntGraphNorm2->SetTitle(";Azimutal Angle #phi [deg.];Reduced Intensity #it{I}_{red}");
    IntGraphNorm2->GetXaxis()->SetLimits(0, 360);
    IntGraphNorm2->GetXaxis()->SetNdivisions(8, 8, 4, false);
    IntGraphNorm2->GetYaxis()->SetRangeUser(0, 45);
    
    //Drawing the Simulation:
    TLine* lowerLine[numberOfPositions];
    TLine* rightLine[numberOfPositions];
    float lineLength;
    for (int i = 0; i < numberOfPositions; i++)
    {
        graphPad1.cd(i + 1);
        gPad->SetLeftMargin(.155);
        gPad->SetBottomMargin(.155);

        //Setting up the lines representing the box walls, such that they are always sufficiently long (same length and both at least reaching the beam position):
        lineLength = 0.22; //For flexible line length depending on position: = Max(MaxElement(1, Gr[i][0]->GetY()), Abs(MinElement(1, Gr[i][0]->GetX())));
        rightLine[i] = new TLine(0, 0, 0, lineLength);
        rightLine[i]->SetLineWidth(5);
        lowerLine[i] = new TLine(-lineLength, 0, 0, 0);
        lowerLine[i]->SetLineWidth(5);
        PosGraph[i]->SetTitle(Form("Position %d;x [m]; y [m]", positionVec[i]));
        GrDummy[i] = new TGraph();
        GrDummy[i]->SetPoint(0, 0.1 * lineLength, -0.1 * lineLength);
        GrDummy[i]->SetPoint(1, -1.1 * lineLength, 1.1 * lineLength);
        GrDummy[i]->SetMarkerColorAlpha(kWhite, 0);
        PosGraph[i]->Add(GrDummy[i]);
        PosGraph[i]->GetXaxis()->SetTitleSize(.055);
        PosGraph[i]->GetXaxis()->SetLabelSize(.045);
        PosGraph[i]->GetYaxis()->SetTitleSize(.055);
        PosGraph[i]->GetYaxis()->SetLabelSize(.045);
        PosGraph[i]->Draw("Y+AP");
        rightLine[i]->Draw("same");
        lowerLine[i]->Draw("same");
        womCircle->Draw("same");
        /*if (i == 1)
        {
            womLabel10->Draw("same");
        }
        else
        {
            womLabel->Draw("same");
        }*/
        womLabel->Draw("same");

        graphPad4.cd(i + 1);
        gPad->SetLeftMargin(.155);
        gPad->SetBottomMargin(.155);
        PosGraph2[i]->SetTitle(Form("Position %d;x [m]; y [m]", positionVec[i]));
        GrDummy2[i] = new TGraph();
        GrDummy2[i]->SetPoint(0, 0.1 * lineLength, -0.1 * lineLength);
        GrDummy2[i]->SetPoint(1, -1.1 * lineLength, 1.1 * lineLength);
        GrDummy2[i]->SetMarkerColorAlpha(kWhite, 0);
        PosGraph2[i]->Add(GrDummy2[i]);
        PosGraph2[i]->GetXaxis()->SetTitleSize(.055);
        PosGraph2[i]->GetXaxis()->SetLabelSize(.045);
        PosGraph2[i]->GetYaxis()->SetTitleSize(.055);
        PosGraph2[i]->GetYaxis()->SetLabelSize(.045);
        PosGraph2[i]->Draw("Y+AP");
        rightLine[i]->Draw("same");
        lowerLine[i]->Draw("same");
        womCircle->Draw("same");
        /*if (i == 1)
        {
            womLabel10->Draw("same");
        }
        else
        {
            womLabel->Draw("same");
        }*/
        womLabel->Draw("same");
    }
    c1->SaveAs(Form("%sSimulation.pdf", out_path.c_str()));
    c4->SaveAs(Form("%sSimulation_noReflection.pdf", out_path.c_str()));
    
    //Drawing the intensities:
    c2->cd();
    c2->SetGrid();
    gPad->SetLeftMargin(.18);
    gPad->SetBottomMargin(.16);
    IntGraph->GetXaxis()->SetTitleSize(.055);
    IntGraph->GetXaxis()->SetLabelSize(.045);
    IntGraph->GetYaxis()->SetTitleSize(.055);
    IntGraph->GetYaxis()->SetLabelSize(.045);
    IntGraph->Draw("apl");
    IntLegend->Draw("same");
    c2->SaveAs(Form("%sIntensities.pdf", out_path.c_str()));

    c3->cd();
    c3->SetGrid();
    gPad->SetLeftMargin(.18);
    gPad->SetBottomMargin(.16);
    IntGraphNorm->GetXaxis()->SetTitleSize(.055);
    IntGraphNorm->GetXaxis()->SetLabelSize(.045);
    IntGraphNorm->GetYaxis()->SetTitleSize(.055);
    IntGraphNorm->GetYaxis()->SetLabelSize(.045);
    IntGraphNorm->Draw("apl");
    IntLegendNorm->Draw("same");
    c3->SaveAs(Form("%sreducedIntensities.pdf", out_path.c_str()));

    c5->cd();
    c5->SetGrid();
    gPad->SetLeftMargin(.18);
    gPad->SetBottomMargin(.16);
    IntGraph2->GetXaxis()->SetTitleSize(.055);
    IntGraph2->GetXaxis()->SetLabelSize(.045);
    IntGraph2->GetYaxis()->SetTitleSize(.055);
    IntGraph2->GetYaxis()->SetLabelSize(.045);
    IntGraph2->Draw("apl");
    IntLegend2->Draw("same");
    c5->SaveAs(Form("%sIntensities_noReflection.pdf", out_path.c_str()));

    c6->cd();
    c6->SetGrid();
    gPad->SetLeftMargin(.18);
    gPad->SetBottomMargin(.16);
    IntGraphNorm2->GetXaxis()->SetTitleSize(.055);
    IntGraphNorm2->GetXaxis()->SetLabelSize(.045);
    IntGraphNorm2->GetYaxis()->SetTitleSize(.055);
    IntGraphNorm2->GetYaxis()->SetLabelSize(.045);
    IntGraphNorm2->Draw("apl");
    IntLegendNorm2->Draw("same");
    c6->SaveAs(Form("%sreducedIntensities_noReflection.pdf", out_path.c_str()));
}