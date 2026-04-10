# ISIS + Claude Code sandbox: the ulimit fix

## TL;DR

**Set `ulimit -n 4096` before running any ISIS camera-loading command from inside Claude Code.** That's it. Every tool works: `campt`, `camrange`, `cam2map`, `caminfo`, `spiceinit`, etc. Output is **100% bit-identical** to the same command run from a normal terminal.

```bash
ulimit -n 4096
campt from=my.cub sample=100 line=100
cam2map from=in.cub to=out.cub map=equi.map
```

## Root cause

Claude Code's sandbox on macOS sets `RLIMIT_NOFILE = INT64_MAX` (effectively unlimited). ISIS 9.0.0's `CubeManager` class sizes its cube cache as a fraction of this value:

```cpp
// isis/src/base/objs/CubeManager/CubeManager.cpp
p_maxOpenFiles = fileLimit.rlim_cur * .60;
```

When `rlim_cur` is `0x7fffffffffffffff`, the multiplication overflows during cast and `p_maxOpenFiles` ends up as a garbage value. This corrupts the internal `QQueue<QString>` that tracks LRU eviction of cached cubes. On the first call to `CubeManager::Open` (triggered from `DemShape::DemShape` when loading the MOLA DEM as a shape model), the corrupted queue's `removeAll` does a `memmove` on an invalid address and the process dies with `SIGSEGV`.

Because ISIS is x86_64 only, the process runs under Rosetta 2 on Apple Silicon. The stack trace reports `_platform_memmove$VARIANT$Rosetta`, which is what initially looked like a sandbox/Rosetta interaction bug but is actually just a consequence of the cache-limit overflow.

## Stack trace (for the record)

```
libsystem_platform.dylib  _platform_memmove$VARIANT$Rosetta
libQt5Core                QListData::remove(int)
libisis9.0.0              Isis::DemShape::DemShape(Isis::Target*, Isis::Pvl&)
libisis9.0.0              Isis::EquatorialCylindricalShape::EquatorialCylindricalShape(...)
libisis9.0.0              Isis::ShapeModelFactory::create(...)
libisis9.0.0              Isis::Target::Target(Isis::Spice*, Isis::Pvl&)
libisis9.0.0              Isis::Spice::init(...)
libisis9.0.0              Isis::Spice::Spice(Isis::Cube&)
libisis9.0.0              Isis::Sensor::Sensor(Isis::Cube&)
libisis9.0.0              Isis::Camera::Camera(Isis::Cube&)
libisis9.0.0              Isis::LineScanCamera::LineScanCamera(Isis::Cube&)
libCTXCamera              Isis::CTXCamera::CTXCamera(Isis::Cube&)
libisis9.0.0              Isis::CameraFactory::Create(Isis::Cube&)
libisis9.0.0              Isis::Cube::camera()
campt main path
```

## Why the user's terminal works fine

A normal macOS terminal starts with `ulimit -n` in the hundreds (256 by default, 1024/4096 after a typical shell startup). `rlim_cur * 0.60` produces a reasonable cache size (e.g. `614`), the cache arithmetic is well-defined, and the QList operations work correctly.

Claude Code raises `RLIMIT_NOFILE` to `INT64_MAX` (possibly to avoid any fd exhaustion issues with the large number of subprocesses / tool invocations in an agent session). That unlimited limit is what breaks ISIS 9.0.0.

## Upstream fix

This is a real ISIS 9.0.0 bug. The correct fix upstream is to clamp `p_maxOpenFiles`:

```cpp
// Clamp to a sane value to avoid overflow when rlim_cur is RLIM_INFINITY.
rlim_t target = (fileLimit.rlim_cur == RLIM_INFINITY)
    ? 1024  // sensible default
    : fileLimit.rlim_cur;
p_maxOpenFiles = std::max<rlim_t>(16, std::min<rlim_t>(target * 0.60, 65536));
```

Worth reporting to DOI-USGS/ISIS3 as a crash under any environment that sets `RLIMIT_NOFILE=RLIM_INFINITY`. This includes Claude Code but also container runtimes with unlimited fd limits.

## Verification

Bit-identical output:

```
Reference (terminal):  27,371,919 valid pixels
Sandbox (ulimit -n 4096): 27,371,919 valid pixels
Intersection:          27,371,919 (100.00%)
Max |diff|:            0.000000
All pixels identical:  27,371,919 / 27,371,919
```

Tested on ISIS 9.0.0 on macOS 15 arm64 with conda-forge build `ha732985_0`, MRO CTX cube `J08_048038_1842_XN_04N287W.lev1.cub`.
