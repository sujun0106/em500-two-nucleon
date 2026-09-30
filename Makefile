# All production Fortran sources live under app/. Build and run at this root.
FC = gfortran
FFLAGS ?= -O2 -fno-fast-math -ffp-contract=off -ffree-line-length-none -frecursive
LDFLAGS ?=
LINALG ?= bundled
LIBS ?= -llapack -lblas
EXE ?= app/fewbody.x
include app/lapack/sources.mk
CORE = app/fb_precision.f90 app/em500_parameters.f90 \
       app/em500/chiral-twobody-potentials.f90 app/em500/chp-set.f90 \
       app/em500_force.f90 app/quadrature.f90 app/deuteron_solver.f90 app/deuteron_observables.f90 \
       app/np_scattering.f90 app/np_observables.f90 \
       app/coulomb/coul90.f app/pp_scattering.f90 app/pp_observables.f90 app/pp_driver.f90
ifeq ($(LINALG),bundled)
MATHSRC = $(LAPACK_SOURCES)
MATHLINK =
else ifeq ($(LINALG),system)
MATHSRC =
MATHLINK = $(LIBS)
else
$(error LINALG must be bundled or system)
endif

.PHONY: all run check check-pp check-deuteron-observables debug clean help verify-sources fetch-external
all: $(EXE)

# The upstream dependency is obtained locally and is excluded from Git.
fetch-external:
	./fetch_external.sh
app/coulomb/coul90.f: | fetch-external
	@test -f "$@"

$(EXE): $(CORE) app/fewbody.f90 $(MATHSRC) Makefile app/lapack/sources.mk
	$(FC) $(FFLAGS) -Japp -Iapp $(LDFLAGS) $(MATHSRC) $(CORE) app/fewbody.f90 $(MATHLINK) -o $@
run: $(EXE)
	./run.sh
check: check-pp
	$(MAKE) check-deuteron-observables
	@test -f app/tests/check_reference.f90 || { echo 'Validation sources missing from app/tests/.'; exit 2; }
	$(FC) $(FFLAGS) -Japp -Iapp $(LDFLAGS) $(MATHSRC) $(CORE) app/tests/check_reference.f90 $(MATHLINK) -o app/check_reference.x
	./app/check_reference.x
	$(FC) $(FFLAGS) -Japp -Iapp $(LDFLAGS) $(MATHSRC) $(CORE) app/tests/check_np_solver.f90 $(MATHLINK) -o app/check_np_solver.x
	./app/check_np_solver.x
	$(FC) $(FFLAGS) -Japp -Iapp $(LDFLAGS) $(MATHSRC) $(CORE) app/tests/check_np_observables.f90 $(MATHLINK) -o app/check_np_observables.x
	./app/check_np_observables.x
verify-sources: app/coulomb/coul90.f
	sha256sum -c app/provenance/em500_sources.sha256
	cd app/lapack && sha256sum -c SHA256SUMS
	cd app/coulomb && sha256sum -c SHA256SUMS
	cd app/tests/reference && sha256sum -c SHA256SUMS

debug:
	$(MAKE) EXE=app/fewbody_debug.x FFLAGS='-O0 -g -fno-fast-math -ffp-contract=off -ffree-line-length-none -frecursive -fcheck=all -fbacktrace' all
clean:
	./remake.sh
help:
	@echo 'make          Fetch/check COUL90, then compile into app/fewbody.x.'
	@echo 'make fetch-external    Download/check the pinned external COUL90 source.'
	@echo './fetch_external.sh --source FILE    Install a verified offline copy.'
	@echo 'make verify-sources   Verify sources; missing COUL90 is downloaded first.'
	@echo './run.sh      Run root input.nml and archive on successful completion.'
	@echo './submit.sh   Submit the same run as a local background process.'
	@echo './run.sh --no-archive   Leave completed runtime files in this root.'
	@echo './archive_run.sh       Archive a completed run manually.'
	@echo 'make check    Run bundled Fortran validation checks.'
	@echo './remove.sh   Delete known root runtime files; preserve inputs and output/.'
	@echo './remake.sh   Remove app/ compilation products (also: make clean).'

# pp adds independent Coulomb/exchange checks; old np/deuteron checks remain.
check-pp: app/coulomb/coul90.f
	$(FC) $(FFLAGS) -Japp -Iapp $(LDFLAGS) $(MATHSRC) $(CORE) app/tests/check_pp_solver.f90 $(MATHLINK) -o app/check_pp_solver.x
	./app/check_pp_solver.x
	$(FC) $(FFLAGS) -Japp -Iapp $(LDFLAGS) $(MATHSRC) $(CORE) app/tests/check_pp_observables.f90 $(MATHLINK) -o app/check_pp_observables.x
	./app/check_pp_observables.x
	@for mode in invalid_zero_angle invalid_pi_angle invalid_eta invalid_q invalid_isospin invalid_pauli invalid_s invalid_nan_angle; do \
	  if ./app/check_pp_observables.x $$mode >/dev/null 2>&1; then echo "Expected rejection: $$mode"; exit 1; fi; \
	done

# Sequential invocation from make check avoids concurrent writes of app/*.mod.
check-deuteron-observables: app/coulomb/coul90.f
	@test -f app/tests/check_deuteron_observables.f90 || { echo 'Validation sources missing from app/tests/.'; exit 2; }
	$(FC) $(FFLAGS) -Japp -Iapp $(LDFLAGS) $(MATHSRC) $(CORE) app/tests/check_deuteron_observables.f90 $(MATHLINK) -o app/check_deuteron_observables.x
	./app/check_deuteron_observables.x
