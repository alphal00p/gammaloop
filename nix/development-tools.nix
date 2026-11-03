{pkgs}:
pkgs.buildEnv {
  name = "gammaloop-development-tools";
  paths = with pkgs; [
    bash
    coreutils
    findutils
    gnutar
    jq
    flock
  ];
}
