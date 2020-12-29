using System.Collections.Generic;
using System.Linq;
using UnityEngine;

namespace DvMod.SteamCutoff
{
    public class FireState
    {
        private static readonly Dictionary<SteamLocoSimulation, FireState> states = new Dictionary<SteamLocoSimulation, FireState>();
        public static FireState Instance(SteamLocoSimulation sim)
        {
            if (states.TryGetValue(sim, out var state))
                return state;
            return states[sim] = new FireState();
        }

        public static FireState? Instance(TrainCar car)
        {
            var sim = car.GetComponent<SteamLocoSimulation>();
            if (sim != null)
                return Instance(sim);
            return null;
        }

        private const float HydrogenAtomicWeight = 1.008f;
        private const float CarbonAtomicWeight = 12.011f;
        private const float OxygenAtomicWeight = 15.999f;
        private const float ExcessOxygenFactor = 1.50f; // needed beyond stoichiometric to ensure full combustion
        private const float CarbonOxygenMassFactor = OxygenAtomicWeight / CarbonAtomicWeight * ExcessOxygenFactor;
        // assume approximately equal molar quantities of carbon and hydrogen
        private const float VolatileOxygenMassFactor =
            (((HydrogenAtomicWeight / (HydrogenAtomicWeight + CarbonAtomicWeight)) / 2f) + // # of oxygen atoms required to burn to H2O
            ((CarbonAtomicWeight / (HydrogenAtomicWeight + CarbonAtomicWeight)) * 2f)) // # of oxygen atoms required to burn to CO2
            * OxygenAtomicWeight * ExcessOxygenFactor;
        private const float CoalDensity = 1346f; // kg/m^3

        /// <summary>Coal chunk radius in m.</summary>
        private const float CoalPieceRadius = 0.02f; // coal passed over 1.25 in = 3.175 cm screen, ~4cm diameter
        private const float CoalPieceVolume = (4f / 3f) * Mathf.PI * CoalPieceRadius * CoalPieceRadius * CoalPieceRadius;
        private const float CoalPieceMass = CoalDensity * CoalPieceVolume;

        /// <summary>Volatile consumption rate (kg/s) per unit surface area (m^2)</summary>
        private const float VolatileConsumptionRateFactor = 0.222f;
        /// <summary>Carbon consumption rate (kg/s) per unit surface area (m^2)</summary>
        private const float CarbonConsumptionRateFactor = VolatileConsumptionRateFactor / 100f;

        private const float CoalChunkMass = 2f; // kg
        private const float PiecesPerChunk = CoalChunkMass / CoalPieceMass;

        private const float CoalCompositionCarbon = 0.6f;
        private const float CoalCompositionVolatile = 0.3f;
        private const float CarbonSpecificEnthalpy = 32.81e3f; // kJ/kg
        private const float VolatileSpecificEnthalpy = CarbonSpecificEnthalpy * 2f; // kJ/kg

        /// <summary>Current oxygen supply as a fraction of oxygen demand.</summary>
        public float oxygenAvailability;
        /// <summary>Mass of carbon in each chunk in kg.</summary>
        public readonly List<float> chunkCarbonMasses = new List<float>();
        /// <summary>Mass of volatiles in each chunk in kg.</summary>
        public readonly List<float> chunkVolatileMasses = new List<float>();

        public void AddCoalChunk()
        {
            chunkCarbonMasses.Insert(0, CoalPieceMass * CoalCompositionCarbon);
            chunkVolatileMasses.Insert(0, CoalPieceMass * CoalCompositionVolatile);
        }

        public float TotalMass() => chunkCarbonMasses.Sum() + chunkVolatileMasses.Sum();

        private const float TwoThirds = 2f / 3f;
        private static readonly float SurfaceAreaFactor = 4f * Mathf.PI * Mathf.Pow(3f / 4f / Mathf.PI, TwoThirds);
        public float PieceSurfaceArea(float coalMass) => SurfaceAreaFactor * Mathf.Pow(coalMass, TwoThirds);

        // oxygen mass required to go from C to CO or CO to CO2
        private float CarbonPartialCombustionRequiredOxygen(float carbonMass) => CarbonOxygenMassFactor * carbonMass;
        private float VolatileCombustionRequiredOxygen(float volatileMass) => VolatileOxygenMassFactor * volatileMass;

        private float PieceCarbonConsumptionRate(float carbonMass) => PieceSurfaceArea(carbonMass) * CarbonConsumptionRateFactor;
        private float PieceVolatileReleaseRate(float carbonMass) => PieceSurfaceArea(carbonMass) * VolatileConsumptionRateFactor;
        private float PieceVolatileConsumptionRate(float carbonMass) => PieceSurfaceArea(carbonMass) * VolatileConsumptionRateFactor;

        private float PieceVolatileReleaseMass(float carbonMass, float volatileMass, float deltaTime) =>
            Mathf.Min(volatileMass, PieceVolatileReleaseRate(carbonMass) * deltaTime);

        // mass of carbon in a piece of mass carbonMass converted to carbon-monoxide
        private float PieceCarbonConsumptionMass(float carbonMass, float oxygenMass, float deltaTime) =>
            Mathf.Min(
                carbonMass,
                Mathf.Min(1f, oxygenMass / CarbonPartialCombustionRequiredOxygen(carbonMass)) *
                    PieceCarbonConsumptionRate(carbonMass) * deltaTime);

        private float VolatileCombustionMass(float volatileMass, float oxygenMass) =>
            volatileMass * Mathf.Min(1f, oxygenMass / VolatileCombustionRequiredOxygen(volatileMass));

        private float CarbonMonoxideCombustionMass(float carbonMass, float oxygenMass) =>
            carbonMass * Mathf.Min(1f, oxygenMass / CarbonPartialCombustionRequiredOxygen(carbonMass));

        /// <summary>Coal surface area in m^3.</summary>
        public float TotalSurfaceArea() => PiecesPerChunk * chunkCarbonMasses.Sum(mass => PieceSurfaceArea(mass / PiecesPerChunk));
        /// <summary>Carbon consumption rate in kg/s with unlimited heat and oxygen.</summary>
        public float MaxCarbonConsumptionRate() => CarbonConsumptionRateFactor * TotalSurfaceArea();
        /// <summary>Volatile consumption rate in kg/s with unlimited heat and oxygen.</summary>
        public float MaxVolatileConsumptionRate() => VolatileConsumptionRateFactor * TotalSurfaceArea();
        /// <summary>Oxygen required for coke combustion to CO.</summary>
        public float COOxygenConsumptionRate() => MaxCarbonConsumptionRate() * CarbonOxygenMassFactor / 2f;
        /// <summary>Maximum oxygen consumption in kg/s.</summary>
        public float MaxOxygenConsumptionRate() =>
            (MaxCarbonConsumptionRate() * CarbonOxygenMassFactor) +
            (MaxVolatileConsumptionRate() * VolatileConsumptionRateFactor);

        public float smoothedHeatYieldRate;
        private float smoothedHeatYieldRateVelo;

        /// <summary>Airflow through stack due to natural convection in kg/s.</summary>
        /// Assuming 2m stack height, 0.5m stack radius, 3m overall height delta, 100 C in smokebox, 20 C at stack outlet.
        /// https://www.engineeringtoolbox.com/natural-draught-ventilation-d_122.html
        private const float PassiveStackFlow = 0.6f;
        /// <summary>Mass ratio of air drawn in vs. high-pressure live or exhaust steam vented.</summary>
        public const float DraftRatio = 2.5f;
        /// <summary>Mass ratio of oxygen in atmospheric air.</summary>
        public const float OxygenRatio = 0.21f;

        /// <summary>Set factors affecting oxygen supply for the fire.</summary>
        /// <param name="exhaustFlow">Amount of steam being exhausted from cylinders and blower in kg/s.</summary>
        /// <returns>Oxygen supply in kg.</returns>
        private float PrimaryOxygenSupply(float exhaustFlow, float damper, float deltaTime) =>
            (PassiveStackFlow + (exhaustFlow * DraftRatio)) * OxygenRatio * damper * deltaTime;
        private float SecondaryOxygenSupply(float exhaustFlow, float fireDoor, float deltaTime) =>
            (PassiveStackFlow + (exhaustFlow * DraftRatio)) * OxygenRatio * fireDoor * deltaTime;

        private const float VolatileIgnitionTemperature = 1000f;
        private (float, float) PieceHeatYield(float carbonMass, float volatileMass, float temperature, float exhaustFlow, float damper, float fireDoor, float deltaTime)
        {
            var oxygenMass = PrimaryOxygenSupply(exhaustFlow, damper, deltaTime);
            var carbonConsumed = PieceCarbonConsumptionMass(carbonMass, oxygenMass, deltaTime);
            var releasedVolatileMass = PieceVolatileReleaseMass(carbonMass, volatileMass, deltaTime);
            var energyReleased = 0.25f * CarbonSpecificEnthalpy * carbonConsumed;
            var oxygenConsumed = CarbonPartialCombustionRequiredOxygen(carbonConsumed);

            oxygenMass += SecondaryOxygenSupply(exhaustFlow, fireDoor, deltaTime);

            if (temperature >= VolatileIgnitionTemperature)
            {
                // generously burn volatiles first
                var combustedVolatileMass = VolatileCombustionMass(releasedVolatileMass, oxygenMass - oxygenConsumed);
                energyReleased += VolatileSpecificEnthalpy * combustedVolatileMass;
                oxygenConsumed += VolatileCombustionRequiredOxygen(combustedVolatileMass);
            }
            // use remaining oxygen to burn CO to CO2
            var carbonMonoxideConsumed = CarbonMonoxideCombustionMass(carbonConsumed, oxygenMass - oxygenConsumed);
            oxygenMass -= CarbonPartialCombustionRequiredOxygen(carbonMonoxideConsumed);
            energyReleased += 0.75f * CarbonSpecificEnthalpy * carbonMonoxideConsumed;

            return (energyReleased, oxygenConsumed / oxygenMass);
        }

        /// <summary>Multiplier on combustion rate based on oxygen availability.</summary>
        public float CombustionMultiplier() => Mathf.Min(oxygenAvailability * 2f, 1f);
        /// <summary>Current coal consumption rate in kg/s.</summary>
        public float CoalConsumptionRate() => CombustionMultiplier() * MaxCoalConsumptionRate();
        /// <summary>Energy yield from coal combustion in kW.</summary>
        public float HeatYieldRate()
        {
            // float co = oxygenAvailability >= 0.5f ? 0.25f : oxygenAvailability / 2f;
            // float co2 = oxygenAvailability >= 0.5f ? (1.5f * oxygenAvailability) - 0.75f : 0f;
            // float combustionEfficiency = co + co2
            float combustionEfficiency = ((oxygenAvailability * oxygenAvailability) + oxygenAvailability) * 0.5f;
            return CoalConsumptionRate() * combustionEfficiency * SpecificEnthalpy * (CoalCompositionCarbon + (2 * CoalCompositionVolatile)) * Main.settings.boilerThermalEfficiency;
        }

        public float SmoothedHeatYieldRate(bool fireOn)
        {
            smoothedHeatYieldRate = Mathf.SmoothDamp(smoothedHeatYieldRate, fireOn ? HeatYieldRate() : 0f, ref smoothedHeatYieldRateVelo, 5f);
            return smoothedHeatYieldRate;
        }

        public (float carbonReleased, float volatilesReleased, float oxygenConsumed) SimulateFirebed(float deltaTime)
        {
            var oxygenRequiredForMaximumRate = chunkCarbonMasses.Sum(mass => CarbonPartialCombustionRequiredOxygen(PieceCarbonConsumptionRate(mass)));
            var carbonReleased = chunkCarbonMasses.Sum(mass => PieceCarbonConsumptionMass())

        }

        public void ConsumeCoal(float deltaTime)
        {
            var radiusChange = deltaTime * CombustionMultiplier() * MaximumRadiusChange;
            for (int i = coalPieceRadii.Count - 1; i >= 0; i--)
            {
                if (coalPieceRadii[i] <= radiusChange)
                    coalPieceRadii.RemoveAt(i);
                else
                    coalPieceRadii[i] -= radiusChange;
            }
        }
    }
}