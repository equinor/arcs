{ pkgs, lib, config, inputs, ... }:

{
  languages.python = {
    enable = true;
    uv.enable = true;
    uv.sync.enable = true;
    venv.enable = true;
  };

  packages = with pkgs; [
    libz
  ];
}
