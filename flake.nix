{
  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/release-24.11"; # that's 23.05
    utils.url = "github:numtide/flake-utils";
    naersk.url = "github:nmattia/naersk";
    naersk.inputs.nixpkgs.follows = "nixpkgs";
    rust-overlay.url = "github:oxalica/rust-overlay";
    rust-overlay.inputs.nixpkgs.follows = "nixpkgs";
  };

  outputs = {
    self,
    nixpkgs,
    utils,
    naersk,
    rust-overlay,
  }:
    utils.lib.eachDefaultSystem (system: let
      #pkgs = nixpkgs.legacyPackages."${system}";
      overlays = [(import rust-overlay)];
      pkgs = import nixpkgs {inherit system overlays;};
      rust = pkgs.rust-bin.stable."1.86.0".default.override {
        targets = ["x86_64-unknown-linux-musl"];
        extensions = ["llvm-tools-preview"];
      };

      # Override the version used in naersk
      naersk-lib = naersk.lib."${system}".override {
        cargo = rust;
        rustc = rust;
      };

      bacon = pkgs.bacon;
    in rec {
      # `nix build`
      packages.mbf-fastq-processor = naersk-lib.buildPackage {
        pname = "mbf_rust_processor";
        root = ./.;
        nativeBuildInputs = with pkgs; [pkg-config];
        buildInputs = with pkgs; [openssl cmake];
        release = true;
        CARGO_PROFILE_RELEASE_debug = "0";
      };
      packages.mbf-fastq-processor_other_linux =
        (naersk-lib.buildPackage {
          pname = "mbf_rust_processor";
          root = ./.;
          nativeBuildInputs = with pkgs; [pkg-config];
          buildInputs = with pkgs; [openssl cmake];
          release = true;
          CARGO_PROFILE_RELEASE_debug = "0";
        })
        .overrideAttrs {
          # make it compatible with other linuxes. It's statically linked anyway
          postInstall = ''
            patchelf $out/bin/mbf_fastq_processor --set-interpreter "/lib64/ld-linux-x86-64.so.2"
          '';
        };
      packages.check = naersk-lib.buildPackage {
        src = ./.;
        mode = "check";
        nativeBuildInputs = with pkgs; [pkg-config zstd.bin];
        buildInputs = with pkgs; [openssl cmake zstd.dev];
      };
      packages.test = naersk-lib.buildPackage {
        src = ./.;
        buildInputs = with pkgs; [openssl cmake];
        mode = "test";
        nativeBuildInputs = with pkgs; [pkg-config cargo-nextest];
        cargoTestCommands = old: ["cargo nextest run $cargo_test_options --no-fail-fast"];
        override = {
          buildPhase = ":";
          postCheck = ''
             # make sure that the friendly panic test outputs a friendly panic
            cargo build --release
             if [ $? -ne 0 ]; then
                 echo "Error: Command failed with non-zero status code"
                 exit 1
             fi

             # Check if stderr contains 'this is embarrasing'
             if grep -q "this is embarrasing" <(echo "$result"); then
                 echo "Error: 'this is embarrasing' found in stderr"
                 exit 1
             fi
          '';
        };
        doCheck = true;
      };
      # haven't been able to get this to work
      # packages.coverage = naersk-lib.buildPackage {
      #   src = ./.;
      #   buildInputs = with pkgs; [openssl cmake];
      #   mode = "test";
      #   nativeBuildInputs = with pkgs; [pkg-config cargo-nextest cargo-llvm-cov];
      #   cargoTestCommands = old: ["cargo llvm-cov nextest --no-tests=fail --run-ignored all"];
      #   override = {
      #     buildPhase = ":";
      #     postCheck = ''
      #       cp  target/llvm-cov/html $out/ -r
      #       '';
      #   };
      #   doCheck = true;
      # };
      #cargoTestCommands = old: ["cargo llvm-cov --html nextest --verbose $cargo_test_options"];

      defaultPackage = packages.mbf-fastq-processor;

      # `nix run`
      apps.mbf-fastq-processor = utils.lib.mkApp {drv = packages.my-project;};
      defaultApp = apps.mbf-fastq-processor;

      # `nix develop`
      devShell = pkgs.mkShell {
        # supply the specific rust version
        nativeBuildInputs = [
          bacon
          pkgs.cargo-audit
          pkgs.cargo-crev
          pkgs.cargo-flamegraph
          pkgs.cargo-insta
          pkgs.cargo-nextest
          pkgs.cargo-llvm-cov
          pkgs.cargo-outdated
          pkgs.cargo-udeps
          pkgs.cargo-vet
          pkgs.cmake
          pkgs.git
          pkgs.openssl
          pkgs.pkg-config
          pkgs.rust-analyzer
          rust
        ];
      };
      devShells.doc = pkgs.mkShell {
        nativeBuildInputs = [pkgs.hugo];
      };
    });
}
# {

