from __future__ import annotations

import sys


def main() -> int:
    try:
        import openmm
        from openmm import Platform, Context, LangevinIntegrator, NonbondedForce, System, Vec3, unit
    except Exception as exc:
        print(f"ERROR: failed to import OpenMM: {exc}", file=sys.stderr)
        return 1

    try:
        import cupy as cp
    except Exception as exc:
        print(f"ERROR: failed to import CuPy: {exc}", file=sys.stderr)
        return 1

    print(f"OpenMM {openmm.__version__}")
    print(f"CuPy {cp.__version__}")

    try:
        device_count = int(cp.cuda.runtime.getDeviceCount())
    except Exception as exc:
        print(f"ERROR: CuPy could not query CUDA devices: {exc}", file=sys.stderr)
        return 1

    if device_count < 1:
        print("ERROR: CuPy reported zero CUDA devices.", file=sys.stderr)
        return 1
    print(f"CuPy devices: {device_count}")

    try:
        cuda_platform = Platform.getPlatformByName("CUDA")
    except Exception as exc:
        print(f"ERROR: OpenMM could not load the CUDA platform: {exc}", file=sys.stderr)
        return 1

    system = System()
    for _ in range(2):
        system.addParticle(39.9)
    force = NonbondedForce()
    force.addParticle(0.0, 0.32, 0.2)
    force.addParticle(0.0, 0.32, 0.2)
    system.addForce(force)

    integrator = LangevinIntegrator(
        300 * unit.kelvin,
        1.0 / unit.picosecond,
        0.002 * unit.picoseconds,
    )
    properties = {"CudaPrecision": "mixed"}

    try:
        context = Context(system, integrator, cuda_platform, properties)
        context.setPositions([Vec3(0.0, 0.0, 0.0), Vec3(0.4, 0.0, 0.0)] * unit.nanometer)
        state = context.getState(getEnergy=True, getForces=True)
        energy = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        first_force = state.getForces(asNumpy=True)[0].value_in_unit(unit.kilojoule_per_mole / unit.nanometer)
        print(f"OpenMM CUDA platform: {cuda_platform.getName()}")
        print(f"Potential energy (kJ/mol): {energy:.6f}")
        print(
            "First force vector (kJ/mol/nm): "
            f"[{first_force[0]:.6f}, {first_force[1]:.6f}, {first_force[2]:.6f}]"
        )
    except Exception as exc:
        print(f"ERROR: OpenMM CUDA context failed: {exc}", file=sys.stderr)
        return 1
    finally:
        try:
            del context
        except Exception:
            pass
        del integrator

    try:
        sanity = float(cp.arange(8, dtype=cp.float32).sum().get())
    except Exception as exc:
        print(f"ERROR: CuPy computation failed: {exc}", file=sys.stderr)
        return 1

    print(f"CuPy sanity sum: {sanity:.1f}")
    print("GPU verification succeeded.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
