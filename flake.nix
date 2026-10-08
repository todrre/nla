{
  description = "NLA mini-project: Python with NumPy, SciPy and Matplotlib";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
  };

  outputs = { self, nixpkgs }:
    let
      system = "x86_64-linux";
      pkgs = nixpkgs.legacyPackages.${system};
      pythonEnv = pkgs.python3.withPackages (ps: [
        ps.numpy
        ps.scipy
        ps.matplotlib
        ps.pillow # image loading (task6.ipynb)
        ps.pyqt6 # GUI backend so plt.show() opens a window
        ps.jupyter # notebooks
      ]);
    in
    {
      # `nix build .#python -o .python-env` gives VS Code a stable interpreter path
      packages.${system}.python = pythonEnv;

      devShells.${system}.default = pkgs.mkShell {
        buildInputs = [
          pythonEnv
          pkgs.qt6.qtwayland # lets the Qt window run natively on Wayland
          pkgs.typst
        ];
      };
    };
}
