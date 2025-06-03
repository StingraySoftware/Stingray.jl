using Test, Stingray
using FITSIO

TEST_DATA_PATH = joinpath(@__DIR__, "_data", "testdata")
@assert isdir(TEST_DATA_PATH)

readdir(TEST_DATA_PATH)

# generate mock data
function mock_data(times, energies; energy_column = "ENERGY")
    test_dir = mktempdir()
    sample_file = joinpath(test_dir, "sample.fits")
    # Create a sample FITS file
    FITS(sample_file, "w") do f
        # Create primary HDU with a small array instead of empty
        # Use a single element array instead of empty
        write(f, [0])
        # Create event table in HDU 2
        table = Dict{String,Array}()
        table["TIME"] = times
        table[energy_column] = energies
        table["INDEX"] = collect(1:length(times))
        write(f, table)
    end
    sample_file
end

mock_times = Float64[1.0, 2.0, 3.0, 4.0, 5.0]
mock_energies = Float64[10.0, 20.0, 15.0, 25.0, 30.0]
sample = mock_data(mock_times, mock_energies)

# try reading mock data
data = readevents(sample)
@test data.times == mock_times
@test data.energies == mock_energies
@test length(data.meta.headers) == 17
@test length(data.meta.extra_columns) == 0

data = readevents(sample; extra_columns = ["INDEX"])
@test length(data.meta.extra_columns) == 1

# try reading sample datasets
data = readevents(joinpath(TEST_DATA_PATH, "monol_testA.evt"))

@test length(data.energies) == 1000
@test length(data.times) == 1000
@test first(data.energies) == 174.0

data = readevents(joinpath(TEST_DATA_PATH, "chandra_test.fits"))
@test length(data.times) == 4612
@test length(data.energies) == 4612
@test data.meta.energy_units == "ENERGY"

# should just read the PI column
data = readevents(joinpath(TEST_DATA_PATH, "chandra_noener_test.fits"))
@test length(data.times) == 4612
@test length(data.energies) == 4612
@test data.meta.energy_units == "PI"

data = readevents(joinpath(TEST_DATA_PATH, "xmm_test.fits"))
@test length(data.times) == 1708244
@test length(data.energies) == 1708244
@test data.meta.energy_units == "PI"

@test_throws AssertionError("Times are not sorted (pass `sort = true` to force sorting).") readevents(
    joinpath(TEST_DATA_PATH, "monol_testA_calib_unsrt.evt"),
)
data = readevents(joinpath(TEST_DATA_PATH, "monol_testA_calib_unsrt.evt"), sort = true)
@test length(data.times) == 1000
@test length(data.energies) == 1000
@test issorted(data.times)
@test data.meta.energy_units == "ENERGY"
# make sure the energy column got sorted correctly
@test data.energies[1:100:1000] ≈
      [9.56, 42.40, 41.36, 3.55, 6.40, 5.55, 40.24, 27.15, 20.31, 35.52] rtol=1e-2
