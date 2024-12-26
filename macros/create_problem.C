void create_problem()
{
    auto file = new TFile("data/data.root");
    auto tree = new TTree("channels");
    tree -> Branch("channel");
    tree -> Branch("energy");
}
