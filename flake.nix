{
  description = "A flake for building a Rust library";

  inputs = {
    nixpkgs.url = "github:nixos/nixpkgs/nixpkgs-unstable";
    flake-utils.url = "github:numtide/flake-utils";
    naersk.url = "github:nix-community/naersk";
    naersk.inputs.nixpkgs.follows = "nixpkgs";
  };

  outputs = { self, nixpkgs, flake-utils, naersk }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        pkgs = nixpkgs.legacyPackages.${system};
        naersk-lib = pkgs.callPackage naersk { };
      in
      {
        devShell = naersk-lib.devShell.override {
          motions = pkgs.lib.attrsets.concatLists [
            naersk-lib.devShell.motions
            [ "cargo" "rustc" "rustfmt" "clippy" ]
          ];
        };
      }
    );
}