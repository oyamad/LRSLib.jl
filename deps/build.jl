using BinDeps

@BinDeps.setup

#lrslib_commit = "2b35cc36c39d1ebca35dbc21a75a35864f0c07b1"
lrslib_commit = "4bb57b5aac78deb1e9a25edf81bea16f5611596c"
#lrsname = "lrslib-061"
lrsname = "lrslib-$lrslib_commit"

# julia installs libgmp10 but not libgmp-dev since it
# does not have to compile C program with GMP,
# e.g. it does not need the headers.
# FIXME BinDeps doesn't work with header, it only works for libraries because it checks for the .so
#       It sees libgmp.so (installed by libgmp10 as a julia dependency) and thinks that it's ok but
#       it does not have the headers
#libgmpdev = library_dependency("libgmp-dev", aliases=["libgmp"])
liblrs = library_dependency("liblrs", aliases=[lrsname, "liblrsgmp"])#, depends=[libgmpdev])

official_repo = "http://cgm.cs.mcgill.ca/~avis/C/lrslib/archive/$lrsname.tar.gz"
#forked_repo = "https://github.com/blegat/lrslib/archive/$lrslib_commit.zip"
forked_repo = "https://github.com/oyamad/lrslib/archive/$lrslib_commit.zip"

#GMP
@linux_only begin
  const has_apt = try success(`apt-get -v`) catch e false end
  const has_yum = try success(`yum --version`) catch e false end
  const has_pacman = try success(`pacman -Qq`) catch e false end
  if has_apt || has_yum
    if has_apt
      pkgname = "libgmp-dev"
      pkgman = "apt-get"
    else
      pkgname = "libgmp-devel or gmp-devel"
      pkgman = "yum"
    end

    println("Warning: The compilation of LRS requires the header gmp.h provided by the package $pkgname.")
    println("If the compilation fails, please install it as follows:")
    println("\$ sudo $pkgman install $pkgname")
  end
end

#LRS
provides(Sources,
        Dict(URI(forked_repo) => liblrs), unpacked_dir="$lrsname")

lrssrcdir = joinpath(BinDeps.srcdir(liblrs), lrsname)
lrsprefixdir = joinpath(BinDeps.usrdir(liblrs))
lrslibdir = joinpath(lrsprefixdir, "lib")

targetdirs = AbstractString["liblrsgmp.$(Libdl.dlext)"]

MAKEFILE = "makefile"
@static if is_apple()
    MAKEFILE = "makefile.osx"
end

provides(BuildProcess,
	(@build_steps begin
		GetSources(liblrs)
		CreateDirectory(lrsprefixdir)
		CreateDirectory(lrslibdir)
		@build_steps begin
			ChangeDirectory(lrssrcdir)
			FileRule(joinpath(lrslibdir,"liblrsgmp.$(Libdl.dlext)"),@build_steps begin
				`make all-shared --file=$MAKEFILE SONAME=liblrsgmp.$(Libdl.dlext).0 SHLINK=liblrsgmp.$(Libdl.dlext)`
				`cp liblrsgmp.$(Libdl.dlext) $lrslibdir/liblrsgmp.$(Libdl.dlext)`
			end)
		end
	end),liblrs)

@BinDeps.install Dict(:liblrs => :liblrs)
