let
  sources = import ./nix/sources.nix;
  pkgs = import sources.nixpkgs { };
  mkShellNewEnv = pkgs.mkShell.override { stdenv = pkgs.gcc10Stdenv; };
  doxygen191 = pkgs.doxygen.overrideAttrs (_: rec {
  name = "doxygen-1.9.1";
  src = pkgs.fetchurl {
    urls = [
      "mirror://sourceforge/doxygen/${name}.src.tar.gz" # faster, with https, etc.
      "http://doxygen.nl/files/${name}.src.tar.gz"
    ];
    sha256 = "1lcif1qi20gf04qyjrx7x367669g17vz2ilgi4cmamp1whdsxbk7";
  };
  });
  eigen339 = pkgs.eigen.overrideAttrs (old: rec {
    version = "3.3.9";
    src = pkgs.fetchFromGitLab {
      owner = "libeigen";
      repo = "eigen";
      rev    = "${version}";
      sha256 = "0m4h9fd5s1pzpncy17r3w0b5a6ywqjajmnr720ndb7fc4bn0dhi4";
    };
    # From https://github.com/foolnotion/aoc2020/blob/master/eigen_include_dir.patch
    patches = [ ./eigen_include_dir.patch ];
  });
in mkShellNewEnv {
  nativeBuildInputs = [ pkgs.cmake ];
  buildInputs = with pkgs; [
    gtest
    bashInteractive
    which
    gdb
    eigen339
    gfortran
    doxygen191
  ];
}
