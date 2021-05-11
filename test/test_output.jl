using Fermi.Error

@testset "Output" begin

    @set printstyle repl
    @set output test_output.dat
    x = @capture_out output("Test Message")

    # Must print message and NOT create file
    @test x == "Test Message\n" && !isfile("test_output.dat")

    @set printstyle both
    x = @capture_out output("Test Message Two!")
    if isfile("test_output.dat")
        y = read("test_output.dat", String)
        # Must print message AND create file
        @test x === "Test Message Two!\n" && y == "Test Message Two!\n" 
        rm("test_output.dat")
    else
        # File was not created properly
        @test false
    end

    @set printstyle file
    x = @capture_out output("Test Message Three!")
    if isfile("test_output.dat")
        y = read("test_output.dat", String)
        # Must print message AND NOT create file
        @test x === "" && y == "Test Message Three!\n" 
        rm("test_output.dat")
    else
        # File was not created properly
        @test false
    end

    @set printstyle none
    x = @capture_out output("Test Message Four!")
    # Must neither print or create a file
    if isfile("test_output.dat") 
        @test false
        rm("test_output.dat")
    else
        @test x == ""
    end

    # Check invalid printstyle
    @set printstyle invalid_printstyle
    @test_throws InvalidFermiOption output("This will never print")

    @reset
end
