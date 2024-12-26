TGraph* fGraphE;
TGraph* fGraphX;
double ProbabilityX(double *x, double *par) { return fGraphX -> Eval(x[0]); }

void run()
{
    gRandom -> SetSeed(time(0));

    auto top = new LKDrawingGroup("top");

    double theta1 = 60;
    double theta2 = 80;
    double theta3 = 100;
    double theta4 = 120; 
    double energyResolution = 0.1;
    double energy_to_adc = 4000/8.;
    const int chMax = 4096;
    const int tbMax = 350;

    //

    auto groupX = top -> CreateGroup();
    fGraphE = (TGraph*) (new TFile("energy_vs_theta.root","read")) -> Get("graph");
    fGraphX = (TGraph*) (new TFile("crosssection_vs_theta.root","read")) -> Get("graph");
    auto probX = new TF1("probX",ProbabilityX,theta1,theta4);
    probX -> SetNpx(1000);
    groupX -> CreateDrawing() -> Add(probX);
    groupX -> CreateDrawing() -> Add(fGraphE);
    groupX -> CreateDrawing() -> Add(fGraphX);

    //

    auto sim = new LKChannelSimulator();
    sim -> SetYMax(chMax);
    sim -> SetTbMax(tbMax);
    sim -> SetNumSmoothing(2);
    sim -> SetSmoothingLength(2);
    sim -> SetPedestalFluctuationLength(4);
    sim -> SetPedestalFluctuationScale(0.10);
    sim -> SetPulseErrorScale(0.05);
    sim -> SetBackGroundLevel(47);

    auto group1 = top -> CreateGroup();
    auto group2 = top -> CreateGroup();

    int numSimulations = 20;
    //for (auto iSim=0; iSim<numSimulations; ++iSim)
    int iSim = 0;
    while (iSim<numSimulations)
    {
        double theta = gRandom -> Uniform(theta1,theta4);
        int detector_id = 0;
        if (theta>theta1&&theta<theta2) detector_id = 1;
        else if (theta>theta3&&theta<theta4) detector_id = 2;
        else continue;
        iSim++;

        double tbHit = 100;
        double energy = fGraphE -> Eval(theta);
        energy = energy + gRandom -> Gaus(energy*energyResolution);
        double amplitude = energy * energy_to_adc;

        auto hist = new TH1D(Form("hist%d",iSim),"",tbMax,0,tbMax);
        hist -> SetMinimum(0);
        hist -> SetStats(0);

        sim -> ResetBuffer();
        sim -> AddFluctuatingPedestal();
        sim -> AddHit(tbHit,amplitude);
        sim -> FillHist(hist);

        LKDrawingGroup *group = nullptr;
        if (detector_id==1) group = group1;
        else if (detector_id==2) group = group2;
        auto draw = group -> CreateDrawing();

        draw -> Add(hist,"","data");
        draw -> AddLegendLine(Form("tbHit=%.2f",tbHit));
        draw -> AddLegendLine(Form("amp=%.2f",amplitude));
        draw -> SetCreateLegend();
    }
    top -> Draw();
}
