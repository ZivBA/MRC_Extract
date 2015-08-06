package meshi.applications.loopBuilding;

public class InternalLoopData {
		protected double evEnergy;
		protected double propEnergy;
		protected double bbHBenergy;
		protected double RMS;

		public InternalLoopData() {};
		
		public InternalLoopData(double evEnergy,
		 double propEnergy,
		 double bbHBenergy,
		 double RMS) {
			this.evEnergy = evEnergy;
			this.propEnergy = propEnergy;
			this.bbHBenergy = bbHBenergy;
			this.RMS = RMS;
		}
		
}
