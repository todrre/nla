{
  description = "NLA mini-project: Python with NumPy, SciPy and Matplotlib";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
  };

  outputs = { self, nixpkgs }:
    let
      system = "x86_64-linux";
      pkgs = nixpkgs.legacyPackages.${system};
    in
    {
      devShells.${system}.default = pkgs.mkShell {
        buildInputs = [
          (pkgs.python3.withPackages (ps: [
            ps.numpy
            ps.scipy
            ps.matplotlib
            ps.pyqt6 # GUI backend so plt.show() opens a window
            ps.jupyter # notebooks (task5.ipynb)
          ]))
          pkgs.qt6.qtwayland # lets the Qt window run natively on Wayland
          pkgs.typst
        ];
      };
    };
}
