void GOSimulation::Run(int numEvents)
{
    for (auto iEvent=0; iEvent<numEvents; ++iEvent)
        RunEvent();
}
