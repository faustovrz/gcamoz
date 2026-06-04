## Genesys auth token bridge
##
## Reuses the SAME browser/Google login that genesysr::user_login() performs,
## but (a) works from non-interactive Rscript runs and (b) persists the token
## (including its refresh token) to disk so later scripts reuse/refresh it
## without opening the browser again.
##
## genesysr gates user_login() behind base `interactive()`, so we replicate its
## flow with httr2 directly. httr2 keys off rlang::is_interactive(), which we
## enable explicitly. The OAuth client is genesysr's built-in public client.

suppressMessages({
  library(genesysr)
  library(httr2)
})

TOKEN_PATH <- path.expand("~/.genesys_token.rds")

.genesys_client <- function() {
  setup_production()
  e <- genesysr:::.genesysEnv
  oauth_client(
    id        = e$client_id,
    secret    = e$client_secret,
    token_url = paste0(e$server, "/oauth/token"),
    auth      = "body"
  )
}

#' Interactive one-time login: opens the browser (Google login), saves token.
genesys_browser_login <- function(redirect_uri = "http://127.0.0.1:48913") {
  options(rlang_interactive = TRUE)   # let httr2 run the browser flow headless-ly
  e <- genesysr:::.genesysEnv
  setup_production()
  token <- oauth_flow_auth_code(
    client       = .genesys_client(),
    auth_url     = paste0(e$server, "/oauth/authorize"),
    scope        = "openid",
    pkce         = TRUE,
    redirect_uri = redirect_uri
  )
  saveRDS(token, TOKEN_PATH)
  message("Saved Genesys token to ", TOKEN_PATH)
  invisible(token)
}

#' Load the saved token, refresh if needed, and arm genesysr for API calls.
genesys_use_saved_token <- function() {
  if (!file.exists(TOKEN_PATH)) {
    stop("No saved token at ", TOKEN_PATH,
         " — run genesys_browser_login() once first.")
  }
  setup_production()
  token <- readRDS(TOKEN_PATH)

  expired <- !is.null(token$expires_at) &&
    as.numeric(token$expires_at) <= as.numeric(Sys.time())

  if (expired && !is.null(token$refresh_token)) {
    message("Access token expired; refreshing via refresh_token ...")
    token <- oauth_flow_refresh(
      client        = .genesys_client(),
      refresh_token = token$refresh_token,
      scope         = "openid"
    )
    if (is.null(token$refresh_token)) token$refresh_token <- readRDS(TOKEN_PATH)$refresh_token
    saveRDS(token, TOKEN_PATH)
  }

  genesysr::authorization(paste("Bearer", token$access_token))
  invisible(token)
}
