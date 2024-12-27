TGraph* fGraphX;
double ProbabilityX(double *x, double *par) { return fGraphX -> Eval(x[0]); }

void run_simulation()
{
    gRandom -> SetSeed(time(0));

    int numSimulations = 100000;
    double energyElectronicsLowLimit = 3.2;
    double energyResolution = 0.02;
    double energyResolutionMeV = 0.3;
    double energy_to_adc = 4000/10.;
    const int chMax = 4096;
    const int tbMax = 350;
    LKBinning bnn_tb(49,321);
    LKBinning bnn_ttc1(60,80);
    LKBinning bnn_ttc2(100,120);
    LKBinning bnn_ttl(160,20,70);
    LKBinning bnn_e(200,0,20);

    auto top = new LKDrawingGroup("top");
    auto groupS = top -> CreateGroup("gS");
    auto group1 = top -> CreateGroup("g1");
    auto group2 = top -> CreateGroup("g2");
    auto groupX = top -> CreateGroup("gX");

    auto graphE = (TGraph*) (new TFile("energy_vs_theta.root","read")) -> Get("graph");
    fGraphX = (TGraph*) (new TFile("crosssection_vs_theta.root","read")) -> Get("graph");
    auto probX = new TF1("probX",ProbabilityX,bnn_ttc1.x1(),bnn_ttc2.x2());
    probX -> SetNpx(1000);
    groupX -> CreateDrawing() -> Add(probX);
    groupX -> CreateDrawing() -> Add(graphE);

    auto chsim = new LKChannelSimulator();
    chsim -> SetYMax(chMax);
    chsim -> SetTbMax(tbMax);
    chsim -> SetNumSmoothing(2);
    chsim -> SetSmoothingLength(2);
    chsim -> SetPedestalFluctuationLength(4);
    chsim -> SetPedestalFluctuationScale(0.10);
    chsim -> SetPulseErrorScale(0.05);
    chsim -> SetBackGroundLevel(47);

    auto histET = (bnn_ttl*bnn_e).NewH2("histET");
    auto histE = bnn_ttl.NewH1("histE");
    groupS -> AddHist(histET);
    groupS -> AddHist(histE);

    int countSim = 0;
    int countElectronicsCut = 0;
    while (countSim<numSimulations)
    {
        int detector_id = 0;
        LKDrawingGroup *group = nullptr;
        double theta = probX -> GetRandom();
        if (bnn_ttc1.IsInside(theta)) { detector_id = 1; group = group1; }
        else if (bnn_ttc2.IsInside(theta)) { detector_id = 2; group = group2; }
        else continue;
        double theta_lab = (180 - theta)/2.;
        if (countSim%10000==0) cout << countSim << endl;
        countSim++;

        double tbHit = bnn_tb.GetRandomUniform();
        double energy = graphE -> Eval(theta);
        energy = energy + gRandom -> Gaus(0,energy*energyResolution);
        energy = energy + gRandom -> Gaus(0,energyResolutionMeV);
        if (energy<energyElectronicsLowLimit) {
            countElectronicsCut++;
            continue;
        }
        double amplitude = energy * energy_to_adc;

        histET -> Fill(theta_lab,energy);
        histE -> Fill(theta_lab);

        chsim -> ResetBuffer();
        chsim -> AddFluctuatingPedestal();
        chsim -> AddHit(tbHit,amplitude);

        if ((detector_id==1 && group1->GetNumDrawings()<16) || (detector_id==2 && group2->GetNumDrawings()<16))
        {
            auto hist = new TH1D(Form("hist%d",countSim),"",tbMax,0,tbMax);
            hist -> SetMinimum(0);
            hist -> SetStats(0);
            chsim -> FillHist(hist);
            auto draw = group -> CreateDrawing(Form("draw%d",countSim));
            draw -> Add(hist,"","data");
            draw -> AddLegendLine(Form("tbHit=%.2f",tbHit));
            draw -> AddLegendLine(Form("amp=%.2f",amplitude));
            draw -> SetCreateLegend();
        }
    }

    top -> Draw("v");
}
