{ pkgs, lib, config, inputs, ... }:

{
  languages.python = {
    enable = true;
    package = (if pkgs.stdenv.isDarwin then pkgs.pkgsx86_64Darwin else pkgs).python312;
    poetry.enable = true;
    poetry.install.enable = true;
  };
}
