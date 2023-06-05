# How to call the new (as of 2023-03-01) ChatGPT API from R
# Get your API key over here: https://platform.openai.com/
#gpt_key <- "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" # Don't share this!

library(httr)
library(stringr)

# Calls the ChatGPT API with the given prompt and returns the answer
ask_chatgpt <- function(prompt) {
  gpt_key <- Sys.getenv("gpt_key")
  response <- POST(
    url = "https://api.openai.com/v1/chat/completions", 
    add_headers(Authorization = paste("Bearer", getOption("gpt_key"))),
    content_type_json(),
    encode = "json",
    body = list(
      model = "gpt-3.5-turbo",
      messages = list(list(
        role = "user", 
        content = prompt
      ))
    )
  )
  str_trim(content(response)$choices[[1]]$message$content)
}