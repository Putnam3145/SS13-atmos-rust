# SS13-atmos-rust
A rust implementation of SS13's gas mixtures. It has a couple public functions, react and share, but their use is not recommended because calling DLLs in byond is extremely slow. Profiling shows reactions to be ~25% as fast when called through a DLL rather than through byond. It's pretty bad.
