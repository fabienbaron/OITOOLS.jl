# Does the QML still parse?
#
# Nothing else in this suite asks. Every other testset here drives the Julia side -- the shell,
# the geometry, the data layers -- and none of it instantiates a QML engine, so a syntax error
# in any .qml file is invisible to all of it. That is not hypothetical: a missing closing brace
# in ObserveTab.qml passed 1195 assertions here and only appeared when the window was opened by
# hand, as `Type ObserveTab unavailable` with the real error one line further down.
#
# The check is the real thing rather than a parser of our own: `gui()` loads Main.qml, which
# instantiates every component file, so anything that stops one of them loading stops this.
# `autoquit_ms` closes the window again a moment later.
#
# IN A SEPARATE PROCESS, deliberately. `src/gui/livecanvas.jl` states the rule the whole GUI is
# built around -- no Makie plot may be created after the QML window exists, because insertion
# allocates GL buffers and the context belongs to Qt's render thread. Opening a window inside
# this process would put every testset that ran after it on the wrong side of that rule, and
# the failure would not look like a QML problem when it came. The cost is a fresh Julia and a
# fresh GLMakie: about 70 s, most of it loading packages rather than loading QML.
@testset "the QML parses, and the window opens" begin
    code = """
    using OITOOLS, GLMakie, QMLMakie, QML
    GUI = Base.get_extension(OITOOLS, :OITOOLSGUIExt)
    gui(GUI.Session(); autoquit_ms = 800)
    """
    out = IOBuffer()
    ok = success(pipeline(`$(Base.julia_cmd()) --project=$(Base.active_project()) -e $code`;
                          stdout = out, stderr = out))
    log = String(take!(out))
    # The engine reports the file and line it gave up on; carry that through rather than
    # failing with a bare `false`, since the point of this test is to say WHERE.
    if !ok
        for line in split(log, '\n')
            occursin(r"\.qml:\d+|Type \w+ unavailable|Failed to load QML", line) &&
                @info "QML load" line
        end
    end
    @test ok
    # A file that loads but whose bindings are broken still warns, and those warnings are the
    # next thing worth catching once this is green.
    @test !occursin("Failed to load QML file", log)
end
